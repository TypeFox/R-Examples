## ------------------------------------------------------------------------
# load the package
library(wrassp)
# get the path to the data that comes with the package
wavPath = system.file('extdata', package='wrassp')
# now list the .wav files so we have some audio files to play with
wavFiles = list.files(wavPath, pattern=glob2rx('*.wav'), full.names=TRUE)

## ------------------------------------------------------------------------
# load an audio file, e.g. the first one in the list above
au = read.AsspDataObj(wavFiles[1])
# show class
class(au)
# print object description
print(au)

## ------------------------------------------------------------------------
# extract duration
dur.AsspDataObj(au)
# extract sampling rate
rate.AsspDataObj(au)
# extract number of records/samples
numRecs.AsspDataObj(au)
# extract additional attributes
attributes(au)

## ------------------------------------------------------------------------
# extract track names
tracks.AsspDataObj(au)
# or an alternative way to extract track names
names(au)
# show head of samples
head(au$audio)

# and we can of course also plot these samples 
# (only plot every 10th element to accelerate plotting)
plot(seq(0,numRecs.AsspDataObj(au) - 1, 10) / rate.AsspDataObj(au), 
     au$audio[c(TRUE, rep(FALSE,9))], 
     type='l', 
     xlab='time (s)', 
     ylab='Audio samples')

## ------------------------------------------------------------------------
# manipulate the audio
au$audio = au$audio * 0.5
# write file to tempdir
dir = tempdir()
writeres = write.AsspDataObj(au, file.path(dir, 'newau.wav'))

## ------------------------------------------------------------------------
# calculate formants and corresponding bandwidth values
fmBwVals = forest(wavFiles[1], toFile=F)
# due to toFile=F this returns an object of the type AsspDataObj and 
# prevents the result being saved to disc as an SSFF file
class(fmBwVals)
# extract track names
# this time the object contains muliple tracks (formants + their bandwidths)
tracks.AsspDataObj(fmBwVals)
# with more than one field (in this case 250 F1/F2/F3/F4 values)
dim(fmBwVals$fm)
# plot the formant values
matplot(seq(0,numRecs.AsspDataObj(fmBwVals) - 1) / rate.AsspDataObj(fmBwVals) + 
          attr(fmBwVals, 'startTime'), 
        fmBwVals$fm, 
        type='l', 
        xlab='time (s)', 
        ylab='Formant frequency (Hz)')

## ------------------------------------------------------------------------
# calculate the fundamental frequency contour
f0vals = ksvF0(wavFiles[1], toFile=F)
# plot the fundamental frequency contour
plot(seq(0,numRecs.AsspDataObj(f0vals) - 1) / rate.AsspDataObj(f0vals) +
       attr(f0vals, 'startTime'),
     f0vals$F0, 
     type='l', 
     xlab='time (s)', 
     ylab='F0 frequency (Hz)')

## ------------------------------------------------------------------------
# calculate the RMS-energy contour for all wavFiles
rmsana(wavFiles, outputDirectory = tempdir())
# list new files using wrasspOutputInfos$rmsana$ext (see below)
rmsFilePaths = list.files(tempdir(), 
                          pattern = paste0('*.',wrasspOutputInfos$rmsana$ext), 
                          full.names = T)
# read first rms file 
rmsvals = read.AsspDataObj(rmsFilePaths[1])
# plot the RMS energy contour
plot(seq(0,numRecs.AsspDataObj(rmsvals) - 1) / rate.AsspDataObj(rmsvals) +
       attr(rmsvals, 'startTime'), 
     rmsvals$rms, 
     type='l', 
     xlab='time (s)', 
     ylab='RMS energy (dB)')

## ------------------------------------------------------------------------
# show all function names
names(wrasspOutputInfos)

## ------------------------------------------------------------------------
# show output infos of function forest
wrasspOutputInfos$forest

