using System;
using System.IO;
using ScottPlot;
using System.Drawing;
using System.Net.Mail;

class SplineInterpolation
{
    static double Func(double x) => x - Math.Sin(x) - 0.25;

    static void Main()
    {
        int[] nodeCounts = { 5, 10, 20, 40 }; // Разные значения n
        int k = 1000; // Количество точек проверки
        double xmin = 0, xmax = 2 * Math.PI;

        string tablePath = "spline_results.txt"; // Файл для сохранения таблицы
        using (StreamWriter writer = new StreamWriter(tablePath))
        {
            writer.WriteLine("Количество узлов (n)\tПроверочные точки (k)\tТип сплайна\tМаксимальное отклонение (RS)");

            foreach (int n in nodeCounts)
            {
                double[] x = new double[n];
                double[] y = new double[n];

                // Заполняем узлы интерполяции
                for (int i = 0; i < n; i++)
                {
                    x[i] = xmin + i * (xmax - xmin) / (n - 1);
                    y[i] = Func(x[i]);
                }

                // Вычисление и построение графиков для каждого типа сплайна
                ProcessSpline(x, y, k, "Линейный", writer);
                ProcessSpline(x, y, k, "Квадратичный", writer);
                ProcessSpline(x, y, k, "Кубический", writer);
            }
        }

        Console.WriteLine($"Результаты сохранены в {tablePath}");
        Console.WriteLine("Файлы сохранены, работа программы завершена, нажмите любую клавишу для выхода");
        Console.ReadKey();
    }

    static void ProcessSpline(double[] x, double[] y, int k, string splineType, StreamWriter writer)
    {
        double xmin = x[0], xmax = x[x.Length - 1];
        double[] xs = new double[k], ysFunc = new double[k], ysSpline = new double[k];
        double maxError = 0;

        // Проверка на k точках
        for (int i = 0; i < k; i++)
        {
            double xi = xmin + i * (xmax - xmin) / (k - 1);
            double fxi = Func(xi);
            double sxi;
            if (splineType == "Линейный")
                sxi = LinearSpline(x, y, xi);
            else if (splineType == "Квадратичный")
                sxi = QuadraticSpline(x, y, xi);
            else if (splineType == "Кубический")
                sxi = CubicSpline(x, y, xi);
            else
                throw new ArgumentException("Неизвестный тип сплайна");
            xs[i] = xi;
            ysFunc[i] = fxi;
            ysSpline[i] = sxi;
            maxError = Math.Max(maxError, Math.Abs(fxi - sxi));
        }

        writer.WriteLine($"{x.Length}\t{k}\t{splineType}\t{maxError:F6}");
        PlotSpline(xs, ysFunc, ysSpline, x.Length, splineType);
    }

    // Линейный сплайн
    static double LinearSpline(double[] x, double[] y, double xVal)
    {
        int i = FindSegment(x, xVal); // Находим интервал
        double t = (xVal - x[i]) / (x[i + 1] - x[i]); //Вычисляем параметр t, который показывыет где находится x искомая на промежутке
        return y[i] * (1 - t) + y[i + 1] * t; //Линейно интерполируем
    }

    // Квадратичный сплайн
    static double QuadraticSpline(double[] x, double[] y, double xVal)
    {
        int i = FindSegment(x, xVal); // Находим интервал
        if (i >= x.Length - 2) i = x.Length - 3; // Ограничиваем индекс, чтобы не выйти за пределы массива

        double h = x[i + 1] - x[i]; 
        double a = y[i]; // Значение в начале интервала
        double b = (y[i + 1] - y[i]) / h; // Приближенный наклон
        double c = (y[i + 2] - 2 * y[i + 1] + y[i]) / (h * h); // учитывает кривизну (вторая производная)
        return a + b * (xVal - x[i]) + c * (xVal - x[i]) * (xVal - x[i + 1]); // Подставляем в уравнение
    }

    // Кубический сплайн (натуральный)
    static double CubicSpline(double[] x, double[] y, double xVal)
    {
        int n = x.Length;
        double[] h = new double[n - 1], alpha = new double[n - 1]; //Длины интервалов h и вспомогательный коэффициент alpha
        double[] l = new double[n], mu = new double[n], z = new double[n]; // Массивы для прогонки
        double[] c = new double[n], b = new double[n - 1], d = new double[n - 1]; // Коэффициенты кубического полинома.
        for (int j = 0; j < n - 1; j++)
            h[j] = x[j + 1] - x[j]; //Вычисление h
        for (int j = 1; j < n - 1; j++)
            alpha[j] = (3 / h[j]) * (y[j + 1] - y[j]) - (3 / h[j - 1]) * (y[j] - y[j - 1]); // Вычисление alpha
        l[0] = 1; mu[0] = 0; z[0] = 0; 
        for (int j = 1; j < n - 1; j++) // Решение системы трехдиагональной матрицы методом прогонки
        {
            l[j] = 2 * (x[j + 1] - x[j - 1]) - h[j - 1] * mu[j - 1]; 
            mu[j] = h[j] / l[j];
            z[j] = (alpha[j] - h[j - 1] * z[j - 1]) / l[j];
        }

        l[n - 1] = 1; z[n - 1] = 0; c[n - 1] = 0; // Обратная прогонка для вычисления коэффициентов
        for (int j = n - 2; j >= 0; j--)
        {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
            d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
        }

        int i = FindSegment(x, xVal); // Вычисление значения y искомого
        double dx = xVal - x[i];
        return y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
    }

    static int FindSegment(double[] x, double xVal)
    {
        int i = 0, j = x.Length - 1;
        while (i < j)
        {
            int mid = (i + j) / 2;
            if (xVal > x[mid])
                i = mid + 1;
            else
                j = mid;
        }
        return Math.Max(0, i - 1);
    }

    static void PlotSpline(double[] xs, double[] ysFunc, double[] ysSpline, int n, string splineType)
    {
        var plt = new ScottPlot.Plot();
        var scatterFunc = plt.Add.Scatter(xs, ysFunc, color: ScottPlot.Colors.Blue);
        var scatterSpline = plt.Add.Scatter(xs, ysSpline, color: ScottPlot.Colors.Red);
        scatterFunc.LegendText = "f(x)";
        scatterSpline.LegendText = $"Сплайн {splineType} (n={n})";
        plt.ShowLegend();
        plt.Title($"Интерполяция: {splineType} (n={n})");
        string fileName = $"spline_{splineType}_n{n}.png".Replace(" ", "_");
        plt.SavePng(fileName, 600, 400);
        Console.WriteLine($"График сохранен: {fileName}");
    }
}
