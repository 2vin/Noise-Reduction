#include "AudioFile.h"	// AudioFile for Audio Handling
#include "AudioFFT.h"	// Audio FFT for Fast Fourier Transform of Audio Signal
#include <math.h>		// log10 

const int bufferSize = 2048;
const int stride = bufferSize/2;
const unsigned int bitRate = 44100;
const int noise_seconds = 3; 

#define PI 3.14159265

AudioFile<float> audioFile;

//  void Example()
//  {
//    const size_t bufferSize = 1024; // Needs to be power of 2!

//    std::vector<float> input(bufferSize, 0.0f);
//    std::vector<float> re(audiofft::AudioFFT::ComplexSize(bufferSize));
//    std::vector<float> im(audiofft::AudioFFT::ComplexSize(bufferSize));
//    std::vector<float> output(bufferSize);

//    audiofft::AudioFFT fft;
//    fft.init(1024);
//    fft.fft(input.data(), re.data(), im.data());
//    fft.ifft(output.data(), re.data(), im.data());
//  }

/* dB Conversion from Amplitude and vice versa */
inline float AmplitudeTodB(float amplitude)
{
  return 20.0f * log10(amplitude);
}

inline float dBToAmplitude(float dB)
{
  return pow(10.0f, dB/20.0f);
}


int main(int argc, char** argv)
{
	audioFile.load (argv[1]);

	int channel = 0;
	int numSamples = audioFile.getNumSamplesPerChannel();

	std::vector<float> dbBuffer(bufferSize);
	
	for (int i = 0; i < numSamples - bufferSize; i=i+stride)
	{
		std::cout<<i/bufferSize<<" ";
		dbBuffer = std::vector<float> (bufferSize, 0.0f);

		for (int j = 0; j < bufferSize; j++)
		{
			float currentSample = audioFile.samples[channel][i+j];
			dbBuffer[j] = AmplitudeTodB(abs(currentSample));
			audioFile.samples[channel][i+j] =  ((currentSample > 0) - (currentSample < 0) )* (dBToAmplitude(dbBuffer[j]- (-dbBuffer[j]*12/60.0)) );
		}
	}

	// float filter = 0.6f, smoothed = audioFile.samples[channel][0];
	// for (int i = 0; i < numSamples; i++)
	// {
	// 	float currentSample = audioFile.samples[channel][i];
	// 	smoothed = filter*currentSample + (1.0 - filter)*smoothed;
	// 	audioFile.samples[channel][i] = smoothed;
	// }

	std::cout<<std::endl;

	// Wave file (implicit)
	audioFile.save (argv[2]);
}