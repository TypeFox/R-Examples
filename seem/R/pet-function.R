pet <- function(eleva, Rad, Temp, RH, wind, plot.out=F){

# psychrometer constant a function of elevation via pressure Kpa
psycro <- 6.8*10^(-3)*(101-0.0115*eleva+5.44*10^(-7)*eleva^2)

# convert to deg Kelvin
Temp.kelvin <- Temp + 273
nT <- length(Temp)

# radiation in MJ/m2
nRad <- length(Rad)
 
PETR <- matrix(nrow=nRad,ncol=nT)
PETR.taylor <- PETR; PETR.penman <- PETR
PETA.penman <- PETR.penman; PET.penman  <- PETR.penman

#priestly-taylor
slope.ptaylor <- exp(21.3 - 5304/(Temp.kelvin))*(5304/(Temp.kelvin^2))
# slope vaporization in Kpa/degC
WR.taylor <- slope.ptaylor/(slope.ptaylor+psycro) 
for(i in 1:nT)
PETR.taylor[,i] <- 30.6*(6.28/365)*Rad*(1-0.23)*WR.taylor[i]

# penman
satura <- 0.1* exp(54.88 - 5.03*log(Temp.kelvin)-6791/Temp.kelvin) #KPa
slope.penman <- (satura/Temp.kelvin)*(6791/Temp.kelvin -5.03)
WR.penman <- slope.penman/(slope.penman+psycro) 
latent <- 2.50 - 0.0022*Temp # MJ/Kg 
for(i in 1:nT) PETR.penman[,i] <- (Rad*(1-0.23)/latent[i])*WR.penman[i]

vap.press <- satura*RH/100
wind.fact <- 2.7 + 1.63*wind 
WA.penman <- 1 - WR.penman
for(i in 1:nT) PETA.penman[,i] <- WA.penman[i]*wind.fact*(satura[i] - vap.press[i])

PET.penman <- PETR.penman + PETA.penman

if (plot.out==T) {
 mat<- matrix(1:4,2,2,byrow=T)
 layout(mat,rep(7/2,2),rep(7/2,2),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

 y <- PETR.penman
 matplot(Rad, y, type="l", lty=1:3, ylab="PET (mm)", xlab="Rad (MJ/m2)",ylim=c(0,10),col=1)
 title("Penman Rad Term",cex.main=0.8);legend("topleft", lty=1:3, legend=paste("T=",Temp,"C"))

 y <- PETA.penman
 matplot(Rad, y, type="l", lty=1:3, ylab="PET (mm)", xlab="Rad (MJ/m2)",ylim=c(0,10),col=1)
 title("Penman Aero Term",cex.main=0.8);legend("topleft", lty=1:3, legend=paste("T=",Temp,"C"))

 y <- PET.penman
 matplot(Rad, y, type="l", lty=1:3, ylab="PET (mm)", xlab="Rad (MJ/m2)",ylim=c(0,10),col=1)
 title("Penman Total",cex.main=0.8); legend("topleft", lty=1:3, legend=paste("T=",Temp,"C"))

 y <- PETR.taylor
 matplot(Rad, y, type="l", lty=1:3, ylab="PET (mm)", xlab="Rad (MJ/m2)",ylim=c(0,10),col=1)
 title("Priestley-Taylor",cex.main=0.8);legend("topleft", lty=1:3, legend=paste("T=",Temp,"C"))

}

out <- data.frame(Rad,PETR.taylor,PET.penman,PETR.penman,PETA.penman)
return(out)

}


