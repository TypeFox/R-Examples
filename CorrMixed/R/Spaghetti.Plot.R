Spaghetti.Plot <- function(Dataset, Outcome, Time, Id, Add.Profiles=TRUE,Add.Mean=TRUE, 
                           Add.Median=FALSE, Col=8, Lwd.Me=3, xlim, ylim, ...){

Object <- x <- Dataset   
Outcome <- x[,paste(substitute(Outcome))]
Time <- x[,paste(substitute(Time))] 
Id <- x[,paste(substitute(Id))]
Data <- data.frame(cbind(Outcome, Time, Id))

max_val_y <- max(Data$Outcome, na.rm=TRUE)
max_val_x <- max(Data$Time, na.rm=TRUE)

if (missing(xlim)==TRUE) {xlim=c(0, max_val_x)}
if (missing(ylim)==TRUE) {ylim=c(0, max_val_y)}

plot(y=Outcome, x=Time, type = "n", ylim=ylim, xlim=xlim, ...)

if (Add.Profiles==TRUE){
for (i in 1: length(unique(Id))){
  lines(y=Data$Outcome[Data$Id==unique(Data$Id)[i]], 
        x=Data$Time[Data$Id==unique(Data$Id)[i]], col=Col)
 }
}

if (Add.Mean==TRUE){
mean <- 
  tapply(Data$Outcome, INDEX = Data$Time, FUN = mean, na.rm=TRUE)
lines(mean, x=unique(Data$Time), lwd=Lwd.Me)
  }

if (Add.Median==TRUE){
  median <- 
    tapply(Data$Outcome, INDEX = Data$Time, FUN = median, na.rm=TRUE)
  lines(median, x=unique(Data$Time), lwd=Lwd.Me)
  }
}