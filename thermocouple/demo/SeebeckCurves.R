# Figure 5-1. Seebeck Curves by Thermocouple Type
# Kerlin, T.W., 1999
# Practical Thermocouple Thermometry
# International Society of Automation (ISA)

xL<-c(-500,2500)
yL<-c(-10,90)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["B"]],type='l',xlim=xL,ylim=yL,main='Seebeck Curves by Thermocouple Type',xlab='Temperature (C)',ylab='Seebeck Coeeficient (uV/C)');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["B"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["B"]])),"B"])[1]
text(n, m,'B');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["E"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["E"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["E"]])),"E"])[1]
text(n, m,'E');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["J"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["J"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["J"]])),"J"])[1]
text(n, m,'J');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["K"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["K"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["K"]])),"K"])[1]
text(n, m,'K');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["N"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["N"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["N"]])),"N"])[1]
text(n, m,'N');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["R"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["R"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["R"]])),"R"])[1]
text(n, m,'R');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["S"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["S"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["S"]])),"S"])[1]
text(n, m,'S');par(new=T)
plot(NominalSeebeckCoefficients[,1], NominalSeebeckCoefficients[["T"]],type='l',xlim=xL,ylim=yL,main='',xlab='',ylab='');par(new=T)
n <- max(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["T"]])),1])
m <- rev(NominalSeebeckCoefficients[which(!is.na(NominalSeebeckCoefficients[["T"]])),"T"])[1]
text(n, m,'T');par(new=T)
