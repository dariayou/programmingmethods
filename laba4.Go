package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"time"
)

func main() {
	rand.Seed(time.Now().UnixNano())

	Fs := 10000
	duration := 1.0
	freqs := [][]float64{
		{100},
		{100, 300, 700},
		{100, 300, 700},
	}
	amps := [][]float64{
		{1},
		{1, 1, 1},
		{3, 2, 1},
	}

	// Задача 1: Генерация сигналов
        fmt.Printf("\n--- Задача 1. генерация спектров---\n")
	for i := 0; i < 3; i++ {
		signal := generateSignal(freqs[i], amps[i], Fs, duration)
		fmt.Printf("\n--- Задание %d ---\nf = %v\na = %v\n", i+1, freqs[i], amps[i])
		fmt.Printf("Первые 5 отсчётов: ")
		for j := 0; j < 5; j++ {
			fmt.Printf("%.2f ", signal[j])
		}
		fmt.Println()
	}

	// Задача 2: Анализ сигналов
        fmt.Printf("\n--- Задача 2. спектр сигнала ---\n")
	for i := 0; i < 3; i++ {
		signal := generateSignal(freqs[i], amps[i], Fs, duration)
		spectrum := computeFFT(signal)

		fmt.Printf("\n--- Задание %d ---\n", i+1)
		fmt.Println("Частота(Гц)\tАмплитуда")
		for k := 0; k < len(spectrum)/4; k++ {
			freq := float64(k) * float64(Fs) / float64(len(signal))
			fmt.Printf("%.1f\t%.4f\n", freq, spectrum[k])
		}

		peaks := findPeaks(spectrum, Fs, 0.1)
		fmt.Printf("\nЗадание %d: основные частоты %.0f Гц\n", i+1, peaks)
		fmt.Printf("Сжатие данных: %d → %d (коэффициент %.1f:1)\n",
			len(signal), len(peaks)*2, float64(len(signal))/float64(len(peaks)*2))
	}

	// Задача 3: Анализ с шумом
	signal := generateSignal(freqs[2], amps[2], Fs, duration)
	fmt.Println("\n--- Задача 3: Влияние шума ---")
	noiseMultipliers := []float64{0.5, 1.0, 2.0, 5.0}
	for _, m := range noiseMultipliers {
		noisySignal := addNoise(signal, m)
		noisySpectrum := computeFFT(noisySignal)
		peaks := findPeaks(noisySpectrum, Fs, 0.1)
		fmt.Printf("Шум ×%.1f: частоты %.0f Гц\n", m, peaks)
	}
}

func generateSignal(freqs, amps []float64, Fs int, duration float64) []float64 {
	N := int(float64(Fs) * duration)
	signal := make([]float64, N)
	for i := 0; i < N; i++ {
		t := float64(i) / float64(Fs)
		sum := 0.0
		for j := range freqs {
			sum += amps[j] * math.Cos(2*math.Pi*freqs[j]*t)
		}
		signal[i] = sum
	}
	return signal
}

func computeFFT(signal []float64) []float64 {
	N := len(signal)
	cSignal := make([]complex128, N)
	for i, v := range signal {
		cSignal[i] = complex(v, 0)
	}
	fftResult := FFT(cSignal)
	spectrum := make([]float64, N)
	for i := 0; i < N; i++ {
		if i == 0 || i == N/2 {
			spectrum[i] = cmplx.Abs(fftResult[i]) / float64(N)
		} else {
			spectrum[i] = 2 * cmplx.Abs(fftResult[i]) / float64(N)
		}
	}
	return spectrum
}

func FFT(x []complex128) []complex128 {
	N := len(x)
	if N <= 1 {
		return x
	}

	even := make([]complex128, N/2)
	odd := make([]complex128, N/2)
	for i := 0; i < N/2; i++ {
		even[i] = x[2*i]
		odd[i] = x[2*i+1]
	}

	even = FFT(even)
	odd = FFT(odd)

	out := make([]complex128, N)
	for k := 0; k < N/2; k++ {
		t := cmplx.Exp(complex(0, -2*math.Pi*float64(k)/float64(N))) * odd[k]
		out[k] = even[k] + t
		out[k+N/2] = even[k] - t
	}
	return out
}

func findPeaks(spectrum []float64, Fs int, threshold float64) []float64 {
	var peaks []float64
	N := len(spectrum)
	for k := 1; k < N/2; k++ {
		freq := float64(k) * float64(Fs) / float64(N)
		if spectrum[k] > threshold {
			peaks = append(peaks, math.Round(freq))
		}
	}
	return peaks
}

func addNoise(signal []float64, multiplier float64) []float64 {
	noisy := make([]float64, len(signal))
	for i := range signal {
		noisy[i] = signal[i] + multiplier*rand.NormFloat64()
	}
	return noisy
}
