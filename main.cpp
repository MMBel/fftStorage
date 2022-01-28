#include <iostream>
#include <fstream>
#include "fftStorage.h"
#include <math.h>

using namespace std;

void    TestFFT(unsigned Frequency, unsigned Amplitude)
{
    unsigned SampleRate = 192000;
    unsigned FrameSize = 32768;
    unsigned Averaging = 1;
    fftStorage fstor;
    fstor.setSampleRate(SampleRate);
    fstor.setFrameSize(FrameSize);
    fstor.setAveraging(Averaging);
    fstor.setFreqLimit(15000);
    fstor.Clear("Gaussian45");
    double  pre2PI = 2.d * M_PI;
    double  SineLength = 1.d*SampleRate/Frequency;
    for(unsigned i=0; i<FrameSize*Averaging; ++i) fstor.AddSample(Amplitude*sin(i*pre2PI/SineLength));
    fstor.DoFFT();
    fpoint MaxP = fstor.GetMaxOptimized(Frequency, fstor.AmountPnt);
    cout << "Found Freq=" << round(MaxP.Frequency) << "\tAmp=" << round(MaxP.Amplitude) << "\tor " << 20*log10(round(MaxP.Amplitude)/32768) << " dB" << endl << endl;
    fstor.SaveStorage("a:\\CPP\\fftStorage\\Storage.txt");
}

void    TestAll(unsigned SampleRate, unsigned FrameSize, unsigned FreqLimit, unsigned Averaging)
{
    struct ShowDev
    {
        ShowDev(fftStorage &fstor, string FName)
        {
            fstor.Clear(FName);
            cout << "Accuracy " << round(fstor.TestDeviation()) << " dB for function " << FName << endl;
        }
    };
    fftStorage fs;
    fs.setSampleRate(SampleRate);
    fs.setAveraging(Averaging);
    fs.setFreqLimit(FreqLimit);
    fs.setFrameSize(FrameSize);
    cout << "Deviation of Amplitude detection:" << endl;
        ShowDev(fs, "Rectangular");
        ShowDev(fs, "Bartlett");
        ShowDev(fs, "Welch");
        ShowDev(fs, "Hamming");
        ShowDev(fs, "Hann");
        ShowDev(fs, "Blackman");
        ShowDev(fs, "BlackmanHarris");
        ShowDev(fs, "Gaussian25");
        ShowDev(fs, "Gaussian35");
        ShowDev(fs, "Gaussian45");
}

int main()
{

TestAll(44100, 32768, 20000, 1);

TestFFT(3000, 16384);



return 0;
}
