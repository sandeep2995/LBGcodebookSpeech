#define framesize 320
#define MaxFrameCount 400
#define Ampscale 5000
#define SpeechLength 600000
#define IgnoreSamples 80
#define samplingrate 16000
#define InitFrames 15
#define inpw "y1.wav" //store the signal in wave form in the file mentioned here
//#define inpt "y1.txt" //store the signal in text form in the file mentioned here
#define duration 3 //time given for recording
#define recmod "C:\\RecordingModule.exe" //path of recording module provided
#define path "C:\\universe\\"
#define p 12
#define precision 6 //precision upto these many number of decimal places
#define delim ' ' //character which you used to sepstral coefficients while storing in file
#define filename "cep.txt" //universe name
#define cbsize 32 // code book size
#define cepsize 12 //number of cepstral coefficients
#define thresholdist 0.05 //threshold distortion in percentage w.r.t previous distortion
