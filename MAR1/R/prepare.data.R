
##  Prepares plankton abundance time-series data for MAR analysis:
	#  1) aggregate by month
	#  2) fill in gaps less than or equal to 'fill.gap' months long via linear interpolation
	#  3) mark blocks of consecutive time-steps
	#  4) change zeros to random number less than 1/2 lowest non-zero value 
		# OR add 1 to all values
	#  5) log10 transform
	#  6) zscore values where z = (value-taxon.mean)/taxon.stdev
		# OR z = value-taxon.increment.mean/taxon.increment.stdev (de-seasons data)



prepare.data <- function(data, increment=c("month","year","week","day"),
					fill.gap=0, replace.0s=c(FALSE,"rand.half","add.ones"), 
					log=FALSE, z.method=c(FALSE,"standard","deseason"),order=NULL) {

##  Argument values:
	#  x = dataframe to be transformed; first column dates, following columns plankton taxa
	#  increment = time-step increment
	#  fill.gap = maximum length of time-step gap to be filled by linear interpolation
	#  replace.0s = "rand.half" OR "add.ones" (in quotes) OR FALSE
			# "rand.half" replaces 0 values with a random value less than 1/2 taxon minimum
			# "add.ones" adds 1 to all values
			# FALSE leaves zeros in data
	#  log = should data be log-transformed?
	#  z.method = "standard" OR "deseason" (in quotes) OR FALSE
			# "standard" subtracts the taxon mean and divides by the taxon standard deviation
			# "deseason" does the same, but the means and standard deviations are month-specific 
					# to remove seasonality
			# FALSE does not apply a z-score to the data
	#  order = date order; "dmy" or "ymd" or "mdy" (must be in quotes); function will attempt
		# to distinguish date order on its own if this argument is not provided

increment<-increment[1]
replace.0s<-replace.0s[1]
z.method<-z.method[1]

# remove any columns that contain only zeros
data1<-data
if(any(apply(data[,-1],2,mean,na.rm=T)==0)){
data1<-data[,-(which(apply(data[,-1],2,mean,na.rm=T)==0)+1)]
cat(paste(length(which(apply(data[,-1],2,mean,na.rm=T)==0)),"columns containing all zeros omitted\n"))
}

#====================================================================================
# 1) FORMAT DATES AND AGGREGATE TAXA BY INCREMENT:
#====================================================================================

dates<-as.character(data1[,1])
d<-strsplit(dates[1],"")[[1]]
sep<-unique(d[is.na(suppressWarnings(as.numeric(d)))])
if(length(sep)!=1) stop("Unrecognized date format. Make sure dates are in first column 
  and use consistent non-numeric seperator between year, month, and day.")
d<-data.frame(apply(t(data.frame(strsplit(dates,sep))),2,as.numeric))

if(is.null(order)){
if(length(which(apply(d,2,max)<13))>1) stop("Cannot distinguish between days and months in dates.
	Provide 'order' argument in function call.")
names(d)[which(apply(d,2,max)>12)]<-"d"
names(d)[which(apply(d,2,max)<13)]<-"m"
names(d)[which(apply(d,2,max)>100)]<-"y"
} else {names(d)<-strsplit(order,"")[[1]]}

data2<-data.frame(date=as.Date(paste(d$y,d$m,d$d,sep="-")),data1[,-1])
data2<-aggregate(data2[,-1],by=list(data2[,1]),mean,na.rm=T)
names(data2)[1]<-"date"
date<-seq.Date(min(data2[,1]),max(data2[,1]),1)
data2<-merge(data.frame(date=date),data2,all=T)
taxa <- data2[,-1]

if(increment=="month") t<-"%Y-%m"
if(increment=="year") t<-"%Y"
if(increment=="week") t<-"%Y-%U"
if(increment=="day") t<-"%Y-%j"

tstep=as.character(date,t)

byinc <- aggregate(data2,by=list(tstep=tstep),mean,na.rm=T)
byinc.taxa<-byinc[,-c(1:2)]

if(increment=="day") incremently<-"daily" else incremently<-paste(increment,"ly",sep="")
cat(paste("data successfully aggregated into",incremently,"time steps\n"))

#====================================================================================
# 2) FIND GAPS AND FILL IN THOSE <= fill.gap:
#====================================================================================

missing<-is.na(apply(byinc.taxa,1,mean))  # where are values missing? 1=yes 0=no

byinc2<-byinc

if(any(missing)&fill.gap>0){

# gets lengths of 1 sequences
lengths.miss <- rle(missing)$lengths[which(rle(missing)$values==T)]
# get indices of all missing values
which.miss <- which(missing)
# get starting row of each missing sequence
start.miss <- which.miss[cumsum(c(1,lengths.miss))[-(length(lengths.miss)+1)]]

# focus on missing sequences less than 'fill.gap' value
start.miss<-start.miss[lengths.miss<=fill.gap]
lengths.miss<-lengths.miss[lengths.miss<=fill.gap]

# replace the 'm' missing value sequences in byinc2 that are <= 'fill.gap' 
	# with linearly interpolated values for each 'j' taxa

if(length(lengths.miss)>0){
for(m in 1:length(lengths.miss)){
for(j in 3:ncol(byinc)){

byinc2[start.miss[m]:(start.miss[m]+lengths.miss[m]-1),j]<-
approx(	x=c(1,lengths.miss[m]+2),	# x=gap length plus 2 endpoints
		# y=abundance values on either side of gap
		y=byinc[c(start.miss[m]-1,start.miss[m]+lengths.miss[m]),j],
		# xout=get interpolation for the missing values
		xout=2:(lengths.miss[m]+1))$y  

}}}

cat(paste(length(lengths.miss),"time gaps <=",fill.gap,"interpolated\n"))

}

#====================================================================================
# 3) MARK BLOCKS OF CONTINUOUS TIME:
#====================================================================================

contin <- 1
for(i in 1:nrow(byinc2)){
if(any(is.na(byinc2[i,])))
contin<-c(contin,contin[i]+1) 
else contin<-c(contin,contin[i])
}

byinc2<-cbind(contin=contin[-1],byinc2) # add contin column to data
if(!any(duplicated(byinc2$contin))){
	cat(paste("NO CONTINUOUS TIME-STEPS IN AGGREGATED SERIES...  Returning",incremently,"data\n\n"))
	return(byinc2[,-2])
	}
byinc3<-na.omit(byinc2)  # strip NA rows from data
cat(paste("Time-series contains",length(which(duplicated(byinc3$contin))),"continuous time-steps\n"))

# remove any columns that contain only zeros
if(any(apply(byinc3[,-c(1:3)],2,mean,na.rm=T)==0))
byinc3<-byinc3[,-(which(apply(byinc3[,-c(1:3)],2,mean,na.rm=T)==0)+3)]

#====================================================================================
# 4) REPLACE ZEROS IN DATA:
#====================================================================================

if(log&replace.0s==F&any(byinc3==0)) 
stop("Zeros must be replaced before data are log-transformed.")

if(replace.0s!=F&any(byinc3<0))
stop("Negative values in data indicate values may already be log-transformed.  Cannot replace zeros.")

no.zeros <- byinc3

#=======================================
	## a)  TO RANDOM NUMBERS LESS THAN 1/2 THE LOWEST VALUE FOR EACH TAXON
#=======================================

if (replace.0s=="rand.half"){

for (j in 4:ncol(no.zeros)){
	nonzeros <- no.zeros[no.zeros[,j]>0,j]
	for (i in 1:nrow(no.zeros)){
		if (no.zeros[i,j]==0) no.zeros[i,j] <- runif(1,0,min(nonzeros)/2) }}
}

#=======================================
	## b)  BY ADDING 1 TO ALL VALUES
#=======================================

if (replace.0s=="add.ones"){
no.zeros <- data.frame(no.zeros[,1:3],no.zeros[,4:ncol(no.zeros)]+1)
}

if(replace.0s!=FALSE) cat(paste(sum(byinc3==0),"zeros replaced\n")) #else 
#cat(paste(sum(byinc3==0),"zeros left in data\n"))

#====================================================================================
# 5) LOG10 TRANSFORM VALUES:
#====================================================================================

if(log&any(no.zeros[,4:ncol(no.zeros)]<0)){
	stop("Data contain negative values. Cannot log-transform.")}

if(log){
logged <- data.frame(no.zeros[,1:3],log10(no.zeros[,4:ncol(no.zeros)]))
cat("data successfully log-transformed\n")} else logged<-no.zeros

#====================================================================================
# 6) Z-SCORE VALUES:
#====================================================================================

#=======================================
	## a)  AND REMOVE SEASONALITY
#=======================================

if(z.method==FALSE) zscores<-logged

if(z.method=="deseason"&increment=="year"){
	warning("Cannot deseason yearly values.  Using standard z-score method.",
			call.=F,immediate.=T)
	z.method<-"standard"}

if (z.method=="deseason"){

# remove year part of tstep to get general increment values
overall.inc<-substr(logged$tstep,6,nchar(logged$tstep))
deseason.zeros<-length(unique(overall.inc))/length(overall.inc)
if(deseason.zeros>.25) warning(paste(
	"Using the deseasoning z-score method has generated a large 
	      number of zeros due to many unique sampling",increment,"values."),
	call.=F,immediate.=F)

# create a data frame of overall increment means 'oim'
oim <- data.frame(aggregate(logged[,-c(1:3)], by=list(increment=overall.inc), mean))

# create a data frame of overall increment standard deviations 'oisd'
oisd <- data.frame(aggregate(logged[,-c(1:3)], by=list(increment=overall.inc), sd))  

oisd[is.na(oisd)]<-1  # sds will be NA if only one occurence of overall value
oisd[oisd==0]<-1  # sds will be 0 if multiple occurences of overall value are equal

zscores <- NULL

for(i in 1:nrow(logged)) {
  inc <- overall.inc[i]
  (logged[i,-c(1:3)]-oim[which(oim$increment==inc),-1])/oisd[which(oisd$increment==inc),-1]->zscore
  zscores<-rbind(zscores,zscore)}

zscores <- data.frame(logged[,1:3],zscores)

cat("data successfully z-scored with deseasoning\n")
}

#=======================================
	## b) WITHOUT REMOVING SEASONALITY
#=======================================

if (z.method=="standard"){

om <- apply(data.frame(logged[,-c(1:3)]),2,mean)  # get overall means for each taxon
osd <- apply(data.frame(logged[,-c(1:3)]),2,sd)  # get overall standard deviations for each taxon

zscores <- logged[,1:3]

for(j in 4:ncol(logged)) {
  zscores<-data.frame(zscores,(logged[,j]-om[j-3])/osd[j-3]) }

names(zscores)<-names(logged)

cat("data successfully z-scored\n")
}

zscores<-zscores[,-2]

zscores
}



