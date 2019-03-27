#include "AudioFile.h"	// AudioFile for Audio Handling
#include "AudioFFT.h"	// Audio FFT for Fast Fourier Transform of Audio Signal
#include <math.h>		// log10 

#define PI 3.14159265

/* Variables declaration */
const float bitRate = 44100;
const int bufferSize = 2048;
const int stride = bufferSize/2;
float sensitivity = 2.0;

int channel = 0;
int numSamples;

const float noise_seconds = 0.13; 
const int noise_samples = float(noise_seconds*bitRate);

AudioFile<float> audioFile;
std::vector<float> sigBuffer, sig_copyBuffer, sig_dbBuffer, sig_hzBufferRe, sig_hzBufferIm;
std::vector<float> noiseBuffer, noise_dbBuffer, noise_hzBufferRe, noise_hzBufferIm;

// Contains Noise Profiles
struct Noise
{
	float mean = 0;
	float std_dev = 0;
};

// Funtion Prototyping
float AmplitudeTodB(float amplitude);
float dBToAmplitude(float dB);

int buffer_AmplitudeTodB(std::vector<float>& amplitude, std::vector<float>& dB)
{
	for(int i=0; i<amplitude.size(); i++)
	{
		dB[i] = AmplitudeTodB(abs(amplitude[i]));
	}
}

int buffer_dBToAmplitude(std::vector<float>& dB, std::vector<float>& amplitude)
{
	for(int i=0; i<dB.size(); i++)
	{
		amplitude[i] = ((amplitude[i]>0)-(amplitude[i]<0))*dBToAmplitude(dB[i]);
	}
}

Noise getNoiseProfile(std::vector<float>& noise_dbBuffer)
{
	Noise noise;
	noise.mean = 0;
	noise.std_dev = 0;

	for(int i=0; i < noise_dbBuffer.size(); i++)
	{
		float abs_dB = abs(noise_dbBuffer[i]);
		if(abs_dB > 100)
			continue;
		noise.mean += abs_dB;
	}
	noise.mean *= 1.0/noise_dbBuffer.size();

	for(int i=0; i < noise_dbBuffer.size(); i++)
	{
		float abs_dB = abs(noise_dbBuffer[i]);
		if(abs_dB > 100)
			continue;
		noise.std_dev += (abs_dB - noise.mean)*(abs_dB - noise.mean);
	}
	noise.std_dev = sqrt(noise.std_dev);
	noise.std_dev *= 1.0/noise_dbBuffer.size();

	
	return noise; 
}

int reduceGain(std::vector<float>& sig_dbBuffer, Noise noise)
{
	for(int i = 0; i<sig_dbBuffer.size(); i++)
	{
		float abs_db = abs(sig_dbBuffer[i]);
		std::cout<<noise.mean+(noise.std_dev*sensitivity)<<std::endl;
		float Gain = 0;
		if(abs_db > noise.mean-(noise.std_dev*sensitivity))
		{
			Gain = -32;
		}
		sig_dbBuffer[i] += Gain; 
	}
}


// int FFT()
// {
//    const size_t fftSize = 1024; // Needs to be power of 2!

//    std::vector<float> input(fftSize, 0.0f);
//    std::vector<float> re(audiofft::AudioFFT::ComplexSize(fftSize));
//    std::vector<float> im(audiofft::AudioFFT::ComplexSize(fftSize));
//    std::vector<float> output(fftSize);

//    audiofft::AudioFFT fft;
//    fft.init(1024);
//    fft.fft(input.data(), re.data(), im.data());
// //    fft.ifft(output.data(), re.data(), im.data());
// }

int main(int argc, char** argv)
{
	// Read WAV file
	audioFile.load (argv[1]);

	numSamples = audioFile.getNumSamplesPerChannel();

	sigBuffer.resize(numSamples-noise_samples);
	sig_copyBuffer.resize(numSamples-noise_samples);
	sig_dbBuffer.resize(numSamples-noise_samples);
	sig_hzBufferRe.resize(numSamples-noise_samples);
	sig_hzBufferIm.resize(numSamples-noise_samples);

	noiseBuffer.resize(noise_samples);
	noise_dbBuffer.resize(noise_samples);
	noise_hzBufferIm.resize(noise_samples);
	noise_hzBufferRe.resize(noise_samples);

	// Get Noise Buffer
	for (int i = 0; i < noise_samples; i++)
	{
		noiseBuffer[i] = (audioFile.samples[channel][i]);
	}
	buffer_AmplitudeTodB(noiseBuffer, noise_dbBuffer);
	// buffer_dBToAmplitude(noise_dbBuffer, noiseBuffer);

	// Get Signal Buffer
	for (int i = noise_samples; i < numSamples; i++)
	{
		sigBuffer[i-noise_samples] = (audioFile.samples[channel][i]);
	}
	sig_copyBuffer = sigBuffer;
	buffer_AmplitudeTodB(sigBuffer, sig_dbBuffer);

	// Get Noise Profile
	Noise nProfile = getNoiseProfile(noise_dbBuffer);

	// Reduce Gain from Noise
	reduceGain(sig_dbBuffer, nProfile);

	buffer_dBToAmplitude(sig_dbBuffer, sigBuffer);

	// Retrieve noise+signal from the buffers
	// for (int i = 0; i < numSamples; i++)
	// {	
	// 	if(i < noise_samples)
	// 	{
	// 		audioFile.samples[channel][i] = noiseBuffer[i];
	// 	}
	// 	else
	// 	{
	// 		audioFile.samples[channel][i] = sigBuffer[i-noise_samples];
	// 	}
	// }
	
	// float filter = 0.6f, smoothed = audioFile.samples[channel][0];
	// for (int i = 0; i < numSamples; i++)
	// {
	// 	float currentSample = audioFile.samples[channel][i];
	// 	smoothed = filter*currentSample + (1.0 - filter)*smoothed;
	// 	audioFile.samples[channel][i] = smoothed;
	// }

	// Wave file (implicit)
	audioFile.save (argv[2]);
}



/* dB Conversion from Amplitude and vice versa */
float AmplitudeTodB(float amplitude)
{
  return 20.0f * log10(amplitude);
}

float dBToAmplitude(float dB)
{
  return pow(10.0f, dB/20.0f);
}
