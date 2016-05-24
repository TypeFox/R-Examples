### R code from vignette source 'RHRV-quickstart.Rnw'

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## install.packages("RHRV", dependencies = TRUE)


###################################################
### code chunk number 2: installDonwnloaded (eval = FALSE)
###################################################
## install.packages("RHRV_XXX",repos=NULL)


###################################################
### code chunk number 3: library
###################################################
library(RHRV)


###################################################
### code chunk number 4: accesingData
###################################################
# HRVData structure containing the heart beats 
data("HRVData")
# HRVData structure storing the results of processing the 
# heart beats: the beats have been filtered, interpolated, ... 
data("HRVProcessedData")


###################################################
### code chunk number 5: creation
###################################################
hrv.data  = CreateHRVData()
hrv.data = SetVerbose(hrv.data, TRUE )


###################################################
### code chunk number 6: RHRV-quickstart.Rnw:157-158
###################################################
options(width=80)


###################################################
### code chunk number 7: loading
###################################################
hrv.data = LoadBeatAscii(hrv.data, "example.beats",
       RecordPath = "beatsFolder")


###################################################
### code chunk number 8: derivating
###################################################
hrv.data = BuildNIHR(hrv.data)


###################################################
### code chunk number 9: filtering
###################################################
hrv.data = FilterNIHR(hrv.data)


###################################################
### code chunk number 10: interpolating
###################################################
# Note that it is not necessary to specify freqhr since it matches with
# the default value: 4 Hz
hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)


###################################################
### code chunk number 11: plottingNIHR
###################################################
PlotNIHR(hrv.data)


###################################################
### code chunk number 12: plottingHR
###################################################
PlotHR(hrv.data)


###################################################
### code chunk number 13: timeAnalysis
###################################################
hrv.data = CreateTimeAnalysis(hrv.data, size = 300,
        interval = 7.8125)


###################################################
### code chunk number 14: completeTimeAnalysis
###################################################
hrv.data = CreateHRVData()
hrv.data = SetVerbose(hrv.data,FALSE)
hrv.data = LoadBeatAscii(hrv.data,"example.beats","beatsFolder")
hrv.data = BuildNIHR(hrv.data)
hrv.data = FilterNIHR(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
hrv.data = CreateTimeAnalysis(hrv.data,size=300,interval = 7.8125)
# We can access "raw" data... let's print separately, the SDNN 
# parameter
cat("The SDNN has a value of ",hrv.data$TimeAnalysis[[1]]$SDNN," msec.\n")


###################################################
### code chunk number 15: creatingFreq
###################################################
hrv.data = CreateFreqAnalysis(hrv.data)


###################################################
### code chunk number 16: STFTanalysis
###################################################
hrv.data = CreateHRVData( )
hrv.data = SetVerbose(hrv.data,FALSE)
hrv.data = LoadBeatAscii(hrv.data,"example.beats","beatsFolder")
hrv.data = BuildNIHR(hrv.data)
hrv.data = FilterNIHR(hrv.data)
hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
# Note that it is not necessary to write the boundaries 
# for the frequency bands, since they match
# the default values
hrv.data = CalculatePowerBand( hrv.data , indexFreqAnalysis= 1,
size = 300, shift = 30, type = "fourier",
ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
LFmin = 0.05, LFmax = 0.15, HFmin = 0.15,   HFmax = 0.4 )


###################################################
### code chunk number 17: STFTanalysis2 (eval = FALSE)
###################################################
## hrv.data = CalculatePowerBand( hrv.data , indexFreqAnalysis= 1,
## size = 300, shift = 30 )


###################################################
### code chunk number 18: waveletAnalysis
###################################################
hrv.data = CreateHRVData( )
hrv.data = SetVerbose(hrv.data,FALSE)
hrv.data = LoadBeatAscii(hrv.data,"example.beats","beatsFolder")
hrv.data = BuildNIHR(hrv.data)
hrv.data = FilterNIHR(hrv.data)
hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
# Note that it is not necessary to write the boundaries
# for the frequency bands, since they match the default values
hrv.data = CalculatePowerBand( hrv.data , indexFreqAnalysis= 1,
 type = "wavelet", wavelet = "la8", bandtolerance = 0.01, relative = FALSE,
ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
 LFmin = 0.05, LFmax = 0.15, HFmin = 0.15,   HFmax = 0.4 )


###################################################
### code chunk number 19: RHRV-quickstart.Rnw:446-452
###################################################
hrv.data = CreateHRVData( )
hrv.data = SetVerbose(hrv.data,FALSE)
hrv.data = LoadBeatAscii(hrv.data,"example.beats","beatsFolder")
hrv.data = BuildNIHR(hrv.data)
hrv.data = FilterNIHR(hrv.data)
hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)


###################################################
### code chunk number 20: bothAnalysis
###################################################
# ...
# create structure, load beats, filter and interpolate
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
# use freqAnalysis number 1 for perfoming 
# Fourier analysis. This time, we do not
# write the band's boundaries
hrv.data = CalculatePowerBand( hrv.data , indexFreqAnalysis= 1,
size = 300, shift = 30, sizesp = 2048, type = "fourier")
# use freqAnalysis number 2 for perfoming 
# wavelet analysis. Note the indexFreqAnalysis = 2!!!
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = CalculatePowerBand( hrv.data , indexFreqAnalysis= 2,
 type = "wavelet", wavelet = "la8", bandtolerance = 0.01, relative = FALSE)


###################################################
### code chunk number 21: plottingFreqFourier
###################################################
# Plotting Fourier analysis
PlotPowerBand(hrv.data, indexFreqAnalysis = 1, ymax = 200, ymaxratio = 1.7)


###################################################
### code chunk number 22: plottingFreqWavelet
###################################################
# Plotting wavelet analysis
PlotPowerBand(hrv.data, indexFreqAnalysis = 2, ymax = 700, ymaxratio = 50)


