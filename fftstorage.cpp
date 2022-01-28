#include <iostream>
#include <fstream>
#include "fftstorage.h"
#include "kiss_fft.c"

fftStorage::~fftStorage()
{
    free(kcfg);
}

void    fftStorage::AddSample(float SampleValue)
{
    IBuffer.push_back(SampleValue);
    if(IBuffer.size()>IBufferSize) IBuffer.pop_front();
}

void    fftStorage::DoFFT()
{
    unsigned CpxSize = mFrameSize;
    kiss_fft_cpx InpBuf[CpxSize];
    kiss_fft_cpx OutBuf[CpxSize];
    for(unsigned av=0; av<mFrameSize/16*mAveraging; av+=mFrameSize/16)
    {
        for(unsigned i=0; i<mFrameSize; i++){ InpBuf[i].r=IBuffer[av+i]; InpBuf[i].i=0; }
        switch (WindowFunc)
        {
            default: break;
        case Bartlett:
        {
            const int nPairs = (mFrameSize - 1)/2;
            const float denom = mFrameSize /2.f;
            InpBuf[0].r = 0.f;
            for(int i=1; i<=nPairs; ++i) { const float value = i/denom; InpBuf[i].r *= value; InpBuf[mFrameSize-i].r *= value; }
        }
        break;
        case Hann:
        {
            const double multiplier = 2 * M_PI / mFrameSize;
            static const double coeff0 = 0.5, coeff1 = -0.5;
            for(unsigned i=0; i<mFrameSize; ++i) InpBuf[i].r *= coeff0 + coeff1*cos(i*multiplier);
        }
        break;
        case Hamming:
        {
            const double multiplier = 2 * M_PI / mFrameSize;
            static const double coeff0 = 0.54, coeff1 = -0.46;
            for(unsigned i=0; i<mFrameSize; ++i) InpBuf[i].r *= coeff0 + coeff1*cos(i*multiplier);
        }
        break;
        case Blackman:
        {
            const double multiplier = 2 * M_PI / mFrameSize;
            const double multiplier2= 2 * multiplier;
            static const double coeff0 = 0.42, coeff1 = -0.5, coeff2 = 0.08;
            for(unsigned i=0; i<mFrameSize; ++i) InpBuf[i].r*= coeff0 + coeff1*cos(i*multiplier) + coeff2*cos(i*multiplier2);
        }
        break;
        case BlackmanHarris:
        {
            const double mult_1 = 2 * M_PI / mFrameSize;
            const double mult_2 = 2 * mult_1;
            const double mult_3 = 3 * mult_1;
            static const double coeff0 = 0.35875, coeff1 = -0.48829, coeff2 = 0.14128, coeff3 = -0.01168;
            for(unsigned i=0; i<mFrameSize; ++i) InpBuf[i].r *= coeff0 + coeff1*cos(i*mult_1) + coeff2*cos(i*mult_2) + coeff3*cos(i*mult_3);
        }
        break;
        case Welch:
        {
            const float N = mFrameSize;
            for(unsigned i=0; i<mFrameSize; ++i) { const float iOverN = i / N; InpBuf[i].r *= 4*iOverN*(1-iOverN); }
        }
        break;
        case Gaussian25:
        {
            static const double A = -2 * 2.5*2.5;
            const float N = mFrameSize;
            for(unsigned i=0; i<mFrameSize; ++i)
            {
                const float iOverN = i/N;
                InpBuf[i].r *= exp(A * (0.25 + (iOverN*iOverN) - iOverN));
            }
        }
        break;
        case Gaussian35:
        {
            static const double A = -2 * 3.5*3.5;
            const float N = mFrameSize;
            for(unsigned i=0; i<mFrameSize; ++i)
            {
                const float iOverN = i/N;
                InpBuf[i].r *= exp(A * (0.25 + (iOverN*iOverN) - iOverN));
            }
        }
        break;
        case Gaussian45:
        {
            static const double A = -2 * 4.5*4.5;
            const float N = mFrameSize;
            for(unsigned i=0; i<mFrameSize; ++i)
            {
                const float iOverN = i/N;
                InpBuf[i].r *= exp(A*(0.25 + (iOverN*iOverN) - iOverN));
            }
        }
        break;
        }
        kiss_fft(kcfg, InpBuf, OutBuf);
        Storage.clear();
        for(unsigned i=0; i<CpxSize; i++)
        {
            float Re=OutBuf[i].r;
            float Im=OutBuf[i].i;
            fpoint pnt;
            pnt.Number=i;
            pnt.Frequency=Resolution*i;
            pnt.Amplitude=float(sqrt(Re*Re+Im*Im));
            if(pnt.Frequency<mFreqLimit) Storage.push_back(pnt);
            else break;
        }
    }
}

void    fftStorage::Clear()
{
    std::string WFunc = "";
    switch (WindowFunc)
    {
        default: break;
        case (Hann):        WFunc="Hann";               break;
        case (Welch):       WFunc="Welch";              break;
        case (Hamming):     WFunc="Hamming";            break;
        case (Blackman):    WFunc="Blackman";           break;
        case (Bartlett):    WFunc="Bartlett";           break;
        case (Gaussian25):  WFunc="Gaussian25";         break;
        case (Gaussian35):  WFunc="Gaussian35";         break;
        case (Gaussian45):  WFunc="Gaussian45";         break;
        case (BlackmanHarris): WFunc="BlackmanHarris";  break;
    }
    if(!Ready) Clear(WFunc);
    else Storage.clear();
}

void    fftStorage::Clear(std::string WindowFuncName)
{
    switch (mSamplRate)
    {
        case (22050): break;
        case (44100): break;
        case (48000): break;
        case (96000): break;
        case (192000): break;
        default: { throw std::runtime_error{"Abnormal Sample Rate!"}; } break;
    }
    switch (mFrameSize)
    {
        case(1024): break;
        case(2048): break;
        case(4096): break;
        case(8192): break;
        case(16384): break;
        case(32768): break;
        default: { throw std::runtime_error{"Abnormal Frame Size!"}; } break;
    }
    if(mFreqLimit==0 || mFreqLimit>mSamplRate/2) mFreqLimit = mSamplRate/2;
    if(mFreqLimit<750) mFreqLimit=750;
    kcfg        = kiss_fft_alloc(mFrameSize, 0, 0, 0);
    Resolution  = 1.d*mSamplRate/mFrameSize;
    mAmountPnt  = mFreqLimit/Resolution;
    IBufferSize = mFrameSize * mAveraging;
    IBuffer.clear();
    Storage.reserve(mAmountPnt*mAveraging+mAveraging);
    WindowFunc=Rectangular;
    if(WindowFuncName=="Hann")            WindowFunc=Hann;
    if(WindowFuncName=="Bartlett")        WindowFunc=Bartlett;
    if(WindowFuncName=="Hamming")         WindowFunc=Hamming;
    if(WindowFuncName=="Blackman")        WindowFunc=Blackman;
    if(WindowFuncName=="BlackmanHarris")  WindowFunc=BlackmanHarris;
    if(WindowFuncName=="Welch")           WindowFunc=Welch;
    if(WindowFuncName=="Gaussian25")      WindowFunc=Gaussian25;
    if(WindowFuncName=="Gaussian35")      WindowFunc=Gaussian35;
    if(WindowFuncName=="Gaussian45")      WindowFunc=Gaussian45;
    unsigned    SineFreq = 750;
    double      pre2PI = 2.d * M_PI;
    double      SineLength = 1.d*mSamplRate/SineFreq;
    unsigned    Amplitude=16384;
    for(size_t i=0; i<mFrameSize*mAveraging; ++i) AddSample(Amplitude*sin(i*pre2PI/SineLength));
    DoFFT();
    fpoint FoundMax = GetPrecizeMax(SineFreq, 1);
    Amplifier = 1.0f/(FoundMax.Amplitude/Amplitude);
    IBuffer.clear();
    Ready=true;
}

fpoint  fftStorage::GetAvgPoint(unsigned Number)
{
    fpoint pnt;
    for(fpoint getp : Storage) if(getp.Number==Number) { pnt.Frequency=getp.Frequency; pnt.Amplitude+=getp.Amplitude; }
    pnt.Amplitude/=mAveraging;
    pnt.Number=Number;
    return pnt;
}

fpoint  fftStorage::GetParabolicPoint(unsigned Number)
{
    using namespace std;
    fpoint outpnt;
    if(Number<2 || Number >mAmountPnt-2) return outpnt;
    fpoint pnt_1 = GetAvgPoint(Number-1);
    fpoint pnt_2 = GetAvgPoint(Number);
    fpoint pnt_3 = GetAvgPoint(Number+1);

    float x1 = pnt_1.Frequency;
    float x2 = pnt_2.Frequency;
    float x3 = pnt_3.Frequency;
    float y1 = pnt_1.Amplitude;
    float y2 = pnt_2.Amplitude;
    float y3 = pnt_3.Amplitude;

    if (y1==0 && y2==0) {
        outpnt.Frequency = x3;
        outpnt.Amplitude = y3;
        return outpnt;
    }

    if (y2==0 && y3==0){
        outpnt.Frequency = x1;
        outpnt.Amplitude = y1;
        return outpnt;
    }

    if (y1==0 && y3==0){
        outpnt.Frequency = x2;
        outpnt.Amplitude = y2;
        return outpnt;
    }
    if (y1>y2) {
        outpnt.Frequency = x1;
        outpnt.Amplitude = y1;
        return outpnt;
    }
    if (y2<y3) {
        outpnt.Frequency = x3;
        outpnt.Amplitude = y3;
        return outpnt;
    }

    float a = (y3-((x3*(y2-y1)+x2*y1-x1*y2)/(x2-x1)))/(x3*(x3-x1-x2)+x1*x2);
    float b = (y2-y1)/(x2-x1)-a*(x1+x2);
    float c = (x2*y1-x1*y2)/(x2-x1)+a*x1*x2;
    float pmax = -b/(a*2);

    outpnt.Frequency = pmax;
    outpnt.Amplitude = a*pmax*pmax + b*pmax + c;

    if(!isnormal(outpnt.Frequency) ||
       !isnormal(outpnt.Amplitude) ||
       abs(outpnt.Amplitude-y2*Amplifier)>500
       ) return pnt_2;

    return outpnt;

}

fpoint  fftStorage::GetMaxNearFreq(unsigned Frequency, unsigned StepsRange)
{
    if (StepsRange==0) StepsRange=1;
    if (StepsRange>mAmountPnt) StepsRange=mAmountPnt;
    fpoint outpnt;
    bool Enable = false;
    float MinFreq = Frequency - Resolution*StepsRange;
    float MaxFreq = Frequency + Resolution*StepsRange;
    for(fpoint pnt : Storage){
        if (pnt.Frequency<=MaxFreq) Enable=true;
        if (pnt.Frequency>MaxFreq || pnt.Frequency<MinFreq) Enable=false;
        if (Enable && outpnt.Amplitude<pnt.Amplitude) outpnt = pnt;
    }
    return outpnt;
}

fpoint  fftStorage::GetPrecizeMax(unsigned Frequency, unsigned StepsRange)
{
    return GetParabolicPoint(GetMaxNearFreq(Frequency, StepsRange).Number);
}

fpoint  fftStorage::GetMaxOptimized(unsigned NearFrequency, unsigned StepsRange)
{
    fpoint outp = GetPrecizeMax(NearFrequency, StepsRange);
    outp.Amplitude*=Amplifier;
    return outp;
}

fpoint  fftStorage::GetPointOptimized(unsigned Number)
{
    fpoint outp = GetAvgPoint(Number);
    outp.Amplitude*=Amplifier;
    return outp;
}

void    fftStorage::setSampleRate(unsigned SampleRate)
{
    mSamplRate=SampleRate;
    Ready=false;
}

void    fftStorage::setFrameSize(unsigned FrameSize)
{
    mFrameSize = FrameSize;
    Ready=false;
}

void    fftStorage::setFreqLimit(unsigned FreqLimit)
{
    mFreqLimit = FreqLimit;
    Ready=false;
}

void    fftStorage::setAveraging(unsigned Averaging)
{
    mAveraging = Averaging;
    Ready=false;
}

void    fftStorage::SaveStorage(std::string Filename)
{
    using namespace std;
    ofstream ofile;
    ofile.open(Filename, ios::trunc);
    if(!ofile.is_open()) { cout << "Not opened" << endl; return; }
    ofile << "Storage Capacity = " << Storage.capacity() << endl;
    ofile << "№\tNum\tFreq\tAmpl" << endl;
    for(unsigned i=0; i<Storage.size(); i++)
    {
        ofile << i
        << "\t" << Storage[i].Number
        << "\t" << round(Storage[i].Frequency)
        << "\t" << round(Storage[i].Amplitude*Amplifier)
        << endl;
    }
    ofile.close();
}

float   fftStorage::TestDeviation()
{
    Clear();
    double  pre2PI = 2.d * M_PI;
    float   Amplitude=16384.f;
    float   Delta = 0.f;
    for(unsigned Freq=10; Freq<mFreqLimit; Freq+=250)
    {
        double  SineLength = 1.d*mSamplRate/Freq;
        for(unsigned i=0; i<mFrameSize*mAveraging; ++i) AddSample(Amplitude*sin((i+SineLength/2)*pre2PI/SineLength));
        DoFFT();
        fpoint pnt = GetMaxOptimized(Freq, 1);
        float nDelta = abs(Amplitude-pnt.Amplitude);
        if (Delta < nDelta && nDelta!=Amplitude) Delta = nDelta;
    }
    Clear();
    return 20*log10(Delta/32767.f);
}
