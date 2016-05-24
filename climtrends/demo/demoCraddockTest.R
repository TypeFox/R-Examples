# Craddock test for testing the homogeneity of yearly average temperatures of Milan vs Turin
# example based on:
# Training session on homogenisation methods (Bologna 17-18 May, 2005)
# Historical Climatology Group
# Istituto di Scienze dell'Atmosfera e del Clima (Institute of Atmospheric Sciences and Climate)
# http://www.isac.cnr.it/climstor/EVENTS/hom_training.html

data(yearly.average.temperature.Turin.Milan)
craddock.result <- CraddockTest(yearly.average.temperature.Turin.Milan,2,3)
plot(yearly.average.temperature.Turin.Milan[,1],craddock.result,type='l',xlab='',ylab='Craddock')
# Example of inhomogeneity
tmp <- yearly.average.temperature.Turin.Milan
tmp[20:43,3] <- tmp[20:43,3]+1 # adding an error of +1 degree to Milan's yearly average temperatures from 1980
craddock.result <- CraddockTest(tmp,2,3)
dev.new()
plot(yearly.average.temperature.Turin.Milan[,1],craddock.result,type='l',xlab='',ylab='Craddock')





