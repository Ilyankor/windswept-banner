%required java file
javaaddpath('C:\Users\miles\Downloads\IRIS-WS-2.0.19.jar');

%fetch the data
source = irisFetch.Traces('AK','PAX','--','BHZ','2021-07-29 6:00:00','2021-07-29 7:00:00');

%data
raw = source.data;

%normalize the amplitude
wave = raw/max(abs(raw));

%write the audio file
audiowrite('earthquake.wav',wave,44100);


