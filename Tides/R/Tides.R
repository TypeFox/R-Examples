TidalCharacteristics <- function (	h,  		#(Water level) time series. data frame with time and h column
					h0 = h$h0, 	#Reference level, either single valued or vector with dimension corresponding to h
					T2 = 5*60*60, 	#'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
					hoffset = 0,   #[Maybe obsolete]Margin on reference level, to cope with small fluctuations in the Water level time series
					filtconst = 1,	#Filtering constant for smoothing the time series
					dtMax = 15,	#maximum accepted time interval in a continuous series. Bigger time intervals are considered to be gaps
					unit = "mins",  # unit of dtMax, Tavg
					Tavg = 12.4*60) #Average period of tidal cycle 
{
if (any(!is.element("time",names(h)))) stop("h must be a data frame with columns time (POSIXt) and h (numeric)")
if (!is.element("POSIXt",class(h$time[1]))) stop("h must be a data frame with columns time (POSIXt) and h (numeric)")
if (is.null(h0)) stop("provide reference level h0")

#Calculate high and low water levels (extrema) in water level time series
output <- extrema(h=h,h0=h0,T2 = T2,filtconst=filtconst, hoffset = hoffset)
HL <- output$HL
tij <- output$h

###
#determine gaps in 'continuous' data
gaps <- gapsts(tij$time,dtMax=dtMax)
if (!is.null(gaps)) 	{gaps$N <- tij$N[match(gaps$t1,tij$time)]	#N counts the tidal cycles
			tij$n <- findInterval(tij$time,gaps$t2)+1	#n counts the continuous series of cycles.	
} else tij$n <- 1

###
#Calculate inundation times and dry times
ITDT <- IT(tij[c("time","h")],h0=tij$h0,dtMax=dtMax, h0marg=hoffset)
ITs <- ITDT$IT
if (!is.null(ITs)) ITs$N <- tij$N[match(ITs$t1,tij$time)]
#When gaps are present in the time series, remove all broken cycles +- 1 cycle number to be sure, i.e. in which a gap of data exists
if (!is.null(gaps)) ITs <- ITs[!is.element(ITs$N, c(gaps$N-1,gaps$N,gaps$N+1)),]
  
DTs <- ITDT$DT
if (!is.null(DTs)) DTs$N <- tij$N[match(DTs$t1,tij$time)]
#When gaps are present in the time series, remove all broken cycles +- 1 cycle number to be sure
if (!is.null(gaps)) DTs <- DTs[!is.element(DTs$N,gaps$N)&!is.element(DTs$N,gaps$N+1)&!is.element(DTs$N,gaps$N - 1),]

###
#Calculate total inundation frequence
if (is.null(gaps)) gapstime <- 0 else gapstime <- sum(unclass(gaps$dt))
Ncycles <- floor(unclass(difftime(max(tij$time,na.rm=T),min(tij$time,na.rm=T),units="mins") - gapstime)/(Tavg))[1]
IF <- IF(HL[HL$HL=="H",]$h,HL[HL$HL=="H",]$h0,N=Ncycles)


TideChars <- list(HL=HL,h=tij,gaps=gaps,IF=IF,ITs=ITs,DTs=DTs,h0 = h0,Ncycles=Ncycles,Tunit=unit)
class(TideChars) <- "Tides" 
return(TideChars)
}

print.Tides <- function(x,...){
  if (is.null(x$DTs$dt)) {cat("Inundation frequency: 100% \n") 
	} else {
	cat("Inundation frequency: ", x$IF, "\n")
	cat("Inundations during time span: ", x$IF*x$Ncycles/100, "\n")}

  if (is.null(x$DTs$dt)) {cat ("WARNING not reliable (IF = 100%). Average inundation height: ", mean(x$HL$h-x$HL$h0),"\n") 
	} else {
	cat("Average inundation height: ", mean(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"]),"\n")
	cat("Average inundation height (per cycle): ", sum(x$HL$h[x$HL$HL=="H"]-x$HL$h0[x$HL$HL=="H"])/x$Ncycles,"\n")}

  if (is.null(x$ITs$dt)) {cat ("Average inundation time: Site never inundated \n") 
	} else {
	cat("Average inundation time: ", mean(x$ITs$dt), x$Tunit, "\n") 
	cat("Average inundation time (per cycle): ", sum(x$ITs$dt)/x$Ncycles, x$Tunit, "\n") 
	cat("Maximal inundation time: ", max(x$ITs$dt), x$Tunit, "\n")
	}


  if (is.null(x$DTs$dt)) {cat ("Average dry time: Site never falls dry \n") 
	} else {  
	cat("Average dry time: ", mean(x$DTs$dt), x$Tunit, "\n")
	cat("Average dry time (per cycle): ", sum(x$DTs$dt)/x$Ncycles, x$Tunit, "\n")
	cat("Maximal dry time: ", max(x$DTs$dt), x$Tunit, "\n")}
  
  cat("Average high water: ", mean(x$HL$h[x$HL$HL=="H"]),"\n")
  cat("Average low water: ", mean(x$HL$h[x$HL$HL=="L"]),"\n")
  
  cat("Time span: ", x$Ncycles, "average (tidal) cycles","\n")

  if (is.null(x$gaps)) {cat("There were no gaps in the time series","\n") 
	} else {
	cat("The time series consists of",max(x$gaps$n),"continuous sub-series","\n")}

}

plot.Tides <- function(x,...){
  plot(x$h$time,x$h$h,type="l",xlab="", ylab="waterlevel",...)
  lines(x$h$time,x$h$h0)
  points(x$HL$time,x$HL$h,col="red",pch=20)
  points(x$HL$time[x$HL$HL=="L"],x$HL$h[x$HL$HL=="L"],col="blue",pch=20)
}


extrema <- function(h,            #(Water level) time series. Data frame with time and h column
                    h0,           #Reference level, either single valued or vector with dimension corresponding to h
                    T2 = 5*60*60, #'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
                    hoffset = 0, #Offset level, to prevent spurious maxima generation due to small fluctuations
                    filtconst = 1 #Filtering constant for smoothing the time series
)
{

h$h0 <- h0

#set all levels < h0 equal to h0
h$ho <- h$h 	#first save original waterlevels
h$h[(h$h)<=h$h0] <- h$h0[(h$h)<=h$h0]

#filter tij data, running average of filtconst (default = 1, no filtering) succesive datapoints,to
#remove small fluctuations

h$hfilt <- filter(h$h,rep(1/filtconst,filtconst))

#remove missing values due to filtering
h <- h[!is.na(h$hfilt),]

#Useful matrix to swith between [H,L] or [T,F] representation of high and low phase of ts
HLTF <- data.frame(HL = c("H","L"), TF = c(TRUE,FALSE))

#Here the core thing happens
#If h[t+T2] < h[t] & h[t-T2] < h[t] then high tide, else low tide
h$TF <- (h$hfilt>approx(x=h$time,y=h$hfilt,xout= pmin(h$time+T2,h$time[length(h$time)]))$y)+hoffset&(h$hfilt>approx(x=h$time,y=h$hfilt,xout=pmax(h$time[1],h$time-T2))$y+hoffset)
h$HL <- HLTF$HL[match(h$TF,HLTF$TF)]

#Give every high and low tide phase a number
h$N <- 0
h$N[2:(length(h$time))] <- 1*(h$HL[1:(length(h$time)-1)] != h$HL[2:(length(h$time))])
h$N[1] <- 1
h$N <- cumsum(h$N)


#Now, find all maxima within each high and low phase
#
#Remark: a very short and clean way of coding would be like this
#
#minmax <- by(h,h$N,function(x,...){
#			switch(x$HL[1], 
#				L = x[which.min(x$h),], 
#				H = x[which.max(x$h),])},
#			simplify=T)
#HL <- do.call(rbind,minmax)
#
#However, this is about twice as slow (due to do.call() and by()) as the following, dirtier code

max <- tapply(h$h,h$N,max)
min <- tapply(h$h,h$N,min)
h$max <- max[h$N]
h$min <- min[h$N]
h$HLval <- 0		
h$HLval <- (h$max==h$h & h$HL=="H")*h$h + (h$min==h$h & h$HL=="L")*h$h   #Pick minimum in low water phase and maximum in high water phase
HL <- h[h$HLval != 0,] 	#Select only High and Low waters
HL <- HL[match(unique(HL$N),HL$N),]	#cleanup: take first occurance in cases where maximum or minimum is reached multiple times during a single water phase.


h$h <- h$ho
return(list(HL = HL[c("time","h","HL","h0")], 		# Data frame with extrema
              h = h[c("time","h","h0","HL","N")]) 	# Original water level data frame with additional attributes
              )
}

gapsts <- function(ts,            	# array of times, consisting of different continuous subseries seperated by large gaps
                    dtMax,        	# maximum time interval in a continuous series, or equivalently minimum interval to be characterized as gap.
                    unit = "mins"  	# unit of dtMax; used when ts belongs to class POSIXt
                    )
{

if (!inherits(ts,"POSIXt")){
timediffs <- ts[1:(length(ts)-1)] - ts[2:(length(ts))]
} else {
timediffs <- difftime(ts[1:(length(ts)-1)], ts[2:(length(ts))],units=unit)
}

if (!any(timediffs < - dtMax)) return(NULL)
#Select gaps > dtMax in a timeseries ts
gaps <- ts[c(timediffs < -dtMax,F)]
gaps <- data.frame(t1 = gaps)
gaps$t2 <- ts[match(gaps$t1,ts)  + 1]
gaps$n <- 1:dim(gaps)[1]

if (!inherits(ts,"POSIXt")){
gaps$dt <- gaps$t2 - gaps$t1
} else {
gaps$dt <- difftime(gaps$t2,gaps$t1,units=unit)
}

return(gaps) #Data frame with the initial time, end time and time difference (unit = unit) of each interval > dtMax
}

IT <- function(h,             #Water level time series. data frame with time and h column
                h0,           #Reference level, either single valued or vector with same length as h
                h0marg = 0.3, #Margin on reference level, to cope with small fluctuations in the Water level time series
                dtMax = 15,   #Maximum time interval in continuous water level series
                unit = "mins"  #Unit of dtMax and output time intervals
                )
{

dry <- h[h$h<=(h0 + h0marg),][c("time","h")]

if (dim(dry)[1] == 0) {		
#If the site never falls dry, inundation time equals the time of the time series
  IT <- data.frame(t1 = h$time[1],t2 = h$time[length(h$time)], dt = difftime(h$time[length(h$time)],h$time[1],units = unit))
  DT <- NULL
} else 
{

wet <- h[h$h>(h0 + h0marg),][c("time","h")] 		#dry time = 'inverse' of inundation time

if (dim(wet)[1] == 0) {	
#If the site is never inundated, dry time equals the time of the time series
  IT <- NULL	
  DT <- data.frame(t1 = h$time[1],t2 = h$time[length(h$time)], dt = difftime(h$time[length(h$time)],h$time[1],units = unit))	
} else
{
if (wet$time[1] > h$time[1]) {
			wet <- rbind(h[1,],wet)
			wet$time[1] <- wet$time[1] - unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (wet$time[length(wet$time)] < h$time[length(h$time)]){
			wet <- rbind(wet,h[length(h$time),])
			wet$time[length(wet$time)] <- wet$time[length(wet$time)] + unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (dry$time[1] > h$time[1]) {
			dry <- rbind(h[1,],dry)
			dry$time[1] <- dry$time[1] - unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (dry$time[length(dry$time)] < h$time[length(h$time)]) {
			dry <- rbind(dry,h[length(h$time),])
			dry$time[length(dry$time)] <- dry$time[length(dry$time)] + unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
# Calculate Inundation Time (IT) as the gaps in the DRY time series
# Calculate Dry Time (DT) as the gaps in the WET time series
# Sum of ITs and DTs should equal total length of time series


IT <- gapsts(dry$time,dtMax,unit=unit)
DT <- gapsts(wet$time,dtMax,unit=unit)

}
}

return(list(IT = IT,     #Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of inundation
            DT = DT))    #Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of dry time
}

IF <- function(H,                  # Vector containing high water levels. 
                h0,                # Reference level for which IF has to be calculated, either single valued or array of length = length(H[1,])
                N  = length(H[,1]) # Number of cycles in time series, equals the number of high water levels when these are complete (= default value)
                )
{
IF <- length(H[H>h0])/N
return(IF*100)               #Inundation frequence [%] at (varying) reference level h0)
}
