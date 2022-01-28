#ifndef FFTSTORAGE_H
#define FFTSTORAGE_H
#include <string>
#include <vector>
#include <deque>
#include "kiss_fft.h"

#ifndef M_PI
#define	M_PI		3.14159265358979323846  /* pi */
#endif


/*
    Mih.M.Belovolov, 2021
    Класс для получения АЧХ посредством БПФ.
    Использование:
        Параметры:
            FreqLimit - действует как фильтр ВЧ, все данные по частотам выше этой отбрасываются.
            FrameSize - размер кадра БПФ.
            Averaging - если 1, то анализируется только 1 кадр. Если 8, то анализируются 8 кадров, сдвинутых
                        каждый относительно предыдущего на 1/16 длины, и к выводу предлагаются точки с
                        усредненными значениями от обработки всех кадров.
            Windowing - Имя оконной функции. Варианты: Rectangular, Hann, Bartlett, Hamming,
                        Blackman, BlackmanHarris, Welch, Gaussian25, Gaussian35, Gaussian45.
                        Применяется к исходным данным перед БПФ для сглаживания боковых лепестков и повышения точности.
                        По умолчанию, или если значение не задано, используется Rectangular.

        Пример использования:
        Создаем экземпляр                   fftStorage fstor;
        Устанавливаем параметры             fstor.setSampleRate(44100);     Частота дискретизации (22050 | 44100 | 48000 | 96000 | 192000)
                                            fstor.setFrameSize(8192);       Размер кадра (1024 | 2048 | 4096 | 8192 | 16384 | 32768)
                                            fstor.setAveraging(4);          Усреднение по ряду кадров
                                            fstor.setFreqLimit(20000);      Лимит частоты
        Калибруем амплификатор и очищаем    fstor.Clear("Hann");
            При этом происходит проверка БПФ на тестовой синусоиде и вычисляется коэффициент амплификации
            для преобразования "магнитуды" в исходную амплитуду.
            Надо заметить, что амплификатор корректно восстановит амплитуду только в случае если сигнал
            занимает весь кадр; а если например только половину его длины, то амплитуда будет уменьшена вдвое.

        Тестируем отклонение точности       fstor.TestDeviation()
            Происходит вычисление амплитуд синусоид ряда частот от 10Гц до FreqLimit и полученные значения
            сравниваются с исходной амплитудой. Функция возвращает максимальное отклонение от нее в децибелах и
            очищает хранилище.
            Величина отклонения зависит от выбранной оконной  функции, а чемпионский титул достается Gaussian45

        Добавляем сэмплы                    fstor.AddSample(FloatValue);
            разумеется, имеет смысл добавлять не больше чем fstor.FrameSize/16*fstor.Averaging+fstor.FrameSize
            сэмплов за раз, в противном случае в приемной очереди начнут стираться входные данные.
            Хотя данные можно просто добавлять в любом количестве и производить БПФ, а старые данные будут
            автоматически удаляться.

        Производим БПФ                      fstor.DoFFT();

        Для вывода используется структура fpoint, содержащая порядковый номер точки, ее частоту и амплитуду.

        Получаем уточненную точку максимума по известной примерной частоте, задав ширину шага в стороны от нее
        Производится уточнение частоты точки максимума путем построения перевернутой параболы к соседним точкам и
        выдача ее вершины в виде fpoint (это означает, что эта точка новая, т.е. в хранилище ее нет)
                                            fstor.GetMaxOptimized(Frequency, 1);
            Примечание. Не всегда существует возможность построить параболу по трем точкам, поэтому в случае
            получения некорректных значений, или если разница между амплитудой ближайшей точки и полученной
            превышает 500 (в диапазоне амплитуд 0-32767), то вместо возможной вершины параболы берется эта точка.

        Получаем все точки последовательно например для построения графика
                                            fpoint pnt;
                                            for(unsigned i=0; i<fstor.AmountPnt; ++i)
                                            {
                                                pnt = fstor.GetPointOptimized(i);
                                                cout << pnt.Frequency << " " pnt.Amplitude << endl;
                                            }

        Сохраняем дамп всех точек           fstor.SaveStorage("c:\\tmp\\Storage.txt");

        Очищаем при необходимости           fstor.Clear();

*/

struct      fpoint{
unsigned    Number=0;
float       Frequency=0,
            Amplitude=0;
};

enum        WinFuncNames : int
{
        Rectangular,
        Hann,
        Bartlett,
        Hamming,
        Blackman,
        BlackmanHarris,
        Welch,
        Gaussian25,
        Gaussian35,
        Gaussian45
};

class fftStorage
{
    unsigned                IBufferSize;
    double                  Resolution;
    float                   Amplifier;
    WinFuncNames            WindowFunc;
    kiss_fft_cfg            kcfg;
    std::vector<fpoint>     Storage;
    std::deque<float>       IBuffer;
    fpoint                  GetParabolicPoint(unsigned Number);
    fpoint                  GetMaxNearFreq(unsigned Frequency, unsigned StepsRange);
    fpoint                  GetPrecizeMax(unsigned Frequency, unsigned StepsRange);
    fpoint                  GetAvgPoint(unsigned Number);
    bool                    Ready = false;
    unsigned                mSamplRate = 44100,
                            mFrameSize = 4096,
                            mFreqLimit = 15000,
                            mAveraging = 1,
                            mAmountPnt = 0;
    std::string             Windowing = "";

    public:
    const unsigned &AmountPnt = mAmountPnt;
    ~fftStorage();
    void    AddSample(float SampleValue);
    void    DoFFT();
    void    Clear();
    void    Clear(std::string WindowFuncName);
    fpoint  GetMaxOptimized(unsigned NearFrequency, unsigned StepsRange);
    fpoint  GetPointOptimized(unsigned Number);
    void    SaveStorage(std::string Filename);
    void    setSampleRate(unsigned SampleRate);
    void    setFrameSize(unsigned FrameSize);
    void    setAveraging(unsigned Averaging);
    void    setFreqLimit(unsigned FrequencyLimit);
    float   TestDeviation();


};


#endif // FFTSTORAGE_H
