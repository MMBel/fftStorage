# fftStorage
C++ класс для получения АЧХ из монофонического звукового файла.
Использует библиотеку <a href="https://github.com/mborgerding/kissfft">kiss-fft</a>

Класс принимает на вход массив сэмплов, производит Быстрое Преобразование Фурье и создает массив с его амплитудно-частотной характеристикой, предоставляя интерфейсы для запросов к этому массиву.

### Порядок использования:
#### Создаем экземпляр:
```fftStorage fstor;```
#### Устанавливаем параметры:
- ```fstor.setSampleRate(44100)```   // Частота дискретизации (22050 | 44100 | 48000 | 96000 | 192000 )
- ```fstor.setFrameSize(8192)```     // Размер кадра БПФ (1024 | 2048 | 4096 | 8192 | 16384 | 32768 )
- ```fstor.setFreqLimit(20000)```    // Фильтр ВЧ, частоты выше 20000 отбрасываются
- ```fstor.setAveraging(4)```        // Преобразовываются 4 кадра со смещением 1/16 относительно предыдущего, каждая запрошенная точка будет являться усредненной

Устанавливаем оконную функцию Ханна, калибруем амлификатор под нее и очищаем приемный буфер:

```fstor.Clear("Hann");```

Происходит проверка БПФ на тестовой синусоиде и вычисляется коэффициент амплификации для преобразования "магнитуды" или амплитудной характеристики БПФ в исходную амплитуду.

 Надо заметить, что амплификатор корректно восстановит амплитуду только если сигнал занимает весь кадр, а если только половину кадра, то полученная амплитуда будет уменьшена вдвое.
 
#### Опционально: 
Тестируем отклонение точности определения исходной амплитуды ```fstor.TestDeviation()```

 Происходит вычисление амплитуд ряда частот от 10 Гц до FreqLimit и полученные значения амплитуд сравниваются с исходной амплитудой.
 
 Функция возвращает максимальное отклонение от нее в децибелах и очищает буфер.
 
 Примечание: величина отклонения зависит от выбранной оконной функции и только от этого.
 
 Стандартный результат:
> - Rectangular -17 dB (оконная функция по-умолчанию)
> - Welch -22 dB 
> - Bartlett -23 dB
> - Hamming -24 dB
> - Hann -25 dB
> - Gaussian25 -25dB
> - Blackman -28 dB
> - BlackmanHarris -30 dB
> - Gaussian35 -30 dB
> - Gaussian45 -34 dB

 Из этого явствует, что наибольшую точность определения амплитуды дает оконная функция <b>Gaussian45</b>
#### Добавляем сэмплы:
```fstor.AddSample(float Value);```

Можно добавлять сэмплы в любом количестве, старые данные из буфера удаляются автоматически.

#### Производим БПФ: 
```fstor.DoFFT();```

После этого шага можем запрашивать данные АЧХ.

Для вывода этих данных используется структура fpoint, содержащая порядковый номер точки в массиве, ее частоту и амплитуду.

#### Получаем предположительную точку максимума по известной примерной частоте, задав ширину шага в стороны от нее (каждый шаг частоты равен SampleRate/FrameSize):

```c
fpoint pnt = fstor.GetMaxOptimized(Frequency, 1);
cout << pnt.Frequency << " " << pnt.Amplitude << endl;
```

Производится уточнение частоты точки максимума путем построения перевернутой параболы к его соседним точкам и выдача вершины в виде fpoint.

Это означает, что эта точка не существует в хранилище, ее существование лишь предполагается. Кроме того, не всегда возможно построить параболу по трем точкам, поэтому в случае получения некорректных значений, или если разница между амплитудой ближайшей точки и этой превышает 500 (в диапазоне амплитуд 0-32767), то вместо возможной вершины параболы берется ближайшая реальная точка.

#### Перебираем все точки последовательно (например для построения графика):
```cpp
fpoint pnt;
for (unsigned  i=0; i<fstor.AmountPnt; ++i)
{
 pnt = fstor.GetPointOptimized(i);
 cout << pnt.Frequency << " " << pnt.Amplitude << endl;
}
```

#### Сохраняем дамп хранилища: 
```fstor.SaveStorage("c:\\Storage.txt");```

#### Очищаем буфер если надо: 
```fstor.Clear();```

Примеры использования см main.cpp
 
 
