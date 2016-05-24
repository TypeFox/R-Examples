#rpath = system.file("extdata",package="spdynmod")
#require(spdynmod)
rpath = paste(find.package('spdynmod'),'/extdata',sep='')

input.Data <- c()
 temp <- read.csv(paste(rpath,'/mc_dynmodv2_7_Data_IPRH.csv',sep=''))
 temp <- list(IPRH = temp)
 input.Data <- c(input.Data, temp)
rm(temp)

inputData <- function(x,name,datalist=input.Data) {
	df=datalist[[name]]
	minT <- min(df[,1],na.rm=T)
	maxT <- max(df[,1],na.rm=T)
	if (x < minT | x > maxT) {
		l <- stats::lm(get(colnames(df)[2])~poly(get(colnames(df)[1]),3),data=df)
		do <- data.frame(x); colnames(do) <- colnames(df)[1]
		o <- predict(l,newdata=do)[[1]]	} else {
	t1 <- max(df[which(df[,1] <= x),1])
	t2 <- min(df[which(df[,1] >= x),1])
	if (t1 == t2) {
		o <- df[t1,2]}
	else {
		w1=1/abs(x-t1);w2=1/abs(x-t2)
	o <- ((df[which( df[,1] == t1),2]*w1)+(df[which( df[,1] == t2),2]*w2)) / (w1+w2) } }
  o }
#----------------
MOD <- function(x,y) {	x %% y }
#----------------
RANDOM <- function(x,y) { runif(1,x,y)}
#----------------
NORMAL <- function(x,y) { rnorm(1,x,y) }
#----------------
POISON <- function(x)  { rpois(1,x) }
#----------------
LOGNORMAL <- function(x,y) { rlnorm(1,x,y) }
#----------------
EXPRAND <- function (x) { rexp(1,x) }
#----------------
SINWAVE <- function(x,y) {
Time<-get('Time')
x * sin(2 * pi * Time / y)
}
#----------------
COSWAVE <- function(x,y) {
Time<-get('Time')
x * cos(2 * pi * Time / y)
}
#----------------
COUNTER <- function(x,y) {
	Time<-get('Time')
	if (Time == time[1]) COUNTER_TEMP <<- x
	if (!exists('COUNTER_TEMP')) COUNTER_TEMP <<- x
	else COUNTER_TEMP <- COUNTER_TEMP  + 1
	if (COUNTER_TEMP == y) COUNTER_TEMP  <<- x
	return(COUNTER_TEMP)}
#--------
TREND <- function(x,y,z=0) {
	DT<-get('DT')
	Time<-get('Time')
	if (!exists('AVERAGE_INPUT')) AVERAGE_INPUT <- z # <<
	CHANGE_IN_AVERAGE <- (x - AVERAGE_INPUT) / y
	AVERAGE_INPUT <- AVERAGE_INPUT + (DT * CHANGE_IN_AVERAGE) # <<
	TREND_IN_INPUT <- (x - AVERAGE_INPUT) / (AVERAGE_INPUT * y)
	if (Time == time[length(time)]) rm(AVERAGE_INPUT,envir=environment(TREND))
	TREND_IN_INPUT}
#-----------------
DELAY <- function(x,y,z=NA) { x } # should be developed!
