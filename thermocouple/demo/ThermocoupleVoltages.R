# Figure 3-1. Thermoelectric EMFs for Standard Thermocouples
# Kerlin, T.W., 1999
# Practical Thermocouple Thermometry
# International Society of Automation (ISA)

xL<-c(-500,2500)
yL<-c(-10,90)
plot(thermocoupleTypeBthermoelectricVoltage[,1], thermocoupleTypeBthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='Thermocouple Voltages',xlab='Temperature (C)',ylab='Thermocouple voltage (mV)');par(new=T)
n <- dim(thermocoupleTypeBthermoelectricVoltage)[1]
text(thermocoupleTypeBthermoelectricVoltage[n,1], thermocoupleTypeBthermoelectricVoltage[n,2],'B');par(new=T)
plot(thermocoupleTypeEthermoelectricVoltage[,1], thermocoupleTypeEthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeEthermoelectricVoltage)[1]
text(thermocoupleTypeEthermoelectricVoltage[n,1], thermocoupleTypeEthermoelectricVoltage[n,2],'E');par(new=T)
plot(thermocoupleTypeJthermoelectricVoltage[,1], thermocoupleTypeJthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeJthermoelectricVoltage)[1]
text(thermocoupleTypeJthermoelectricVoltage[n,1], thermocoupleTypeJthermoelectricVoltage[n,2],'J');par(new=T)
plot(thermocoupleTypeKthermoelectricVoltage[,1], thermocoupleTypeKthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeKthermoelectricVoltage)[1]
text(thermocoupleTypeKthermoelectricVoltage[n,1], thermocoupleTypeKthermoelectricVoltage[n,2],'K');par(new=T)
plot(thermocoupleTypeNthermoelectricVoltage[,1], thermocoupleTypeNthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeNthermoelectricVoltage)[1]
text(thermocoupleTypeNthermoelectricVoltage[n,1], thermocoupleTypeNthermoelectricVoltage[n,2],'N');par(new=T)
plot(thermocoupleTypeRthermoelectricVoltage[,1], thermocoupleTypeRthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeRthermoelectricVoltage)[1]
text(thermocoupleTypeRthermoelectricVoltage[n,1], thermocoupleTypeRthermoelectricVoltage[n,2],'R');par(new=T)
plot(thermocoupleTypeSthermoelectricVoltage[,1], thermocoupleTypeSthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeSthermoelectricVoltage)[1]
text(thermocoupleTypeSthermoelectricVoltage[n,1], thermocoupleTypeSthermoelectricVoltage[n,2],'S');par(new=T)
plot(thermocoupleTypeTthermoelectricVoltage[,1], thermocoupleTypeTthermoelectricVoltage[,2],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- dim(thermocoupleTypeTthermoelectricVoltage)[1]
text(thermocoupleTypeTthermoelectricVoltage[n,1], thermocoupleTypeTthermoelectricVoltage[n,2],'T');par(new=T)

