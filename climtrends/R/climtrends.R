



MonthlyAnomaliesFromDailyData<-function(dataSeries,yearAnomalies=NA,fromYear=1961,toYear=1990)
{
if (is.na(yearAnomalies)) stop('<<<yearAnomalies>>> must be a valid year')
monthClimat<-MonthlyClimatologyFromDailyData(dataSeries)
tmpdaily1year<-ValuesBetween2years(dataSeries,yearAnomalies,yearAnomalies)
colnames(tmpdaily1year)<-c('date','data')
monthMean1year<-MonthlyFuncFromDay(tmpdaily1year)
monthMean1year<-monthMean1year[,-1]
monthlyAnomalies<-cbind(month=1:12,NA)
monthlyAnomalies[,2]<-monthMean1year[,2]-MonthlyClimatologyFromDailyData(dataSeries,fromYear,toYear)[,2]
monthlyAnomalies
}

MonthlyClimatologyFromDailyData<-function(dataSeries,fromYear=1961,toYear=1990)
{# monthly climatology (long-term average, for each month, of a given variable)
var61.90<-ValuesBetween2years(dataSeries,fromYear,toYear)
mvar61.90<-MonthlyFuncFromDay(var61.90)
mean.month61.90<-aggregate(data ~ month, mvar61.90, mean)
mean.month61.90
}

DailyClimatologyFromDailyData<-function(dataSeries,fromYear=1961,toYear=1990)
{ # daily climatology (long-term average, for each day, of a given variable)
var61.90<-ValuesBetween2years(dataSeries,fromYear,toYear)
colnames(var61.90)<-c('date','data')
mean.day61.90<-aggregate(data ~ strptime(date, format="%Y-%m-%d")$yday + 1, var61.90, mean)
mean.day61.90
}

ListContinuousDays<-function(dataSeries)
{#show continuous days in daily data
dataSeries[which(!is.na(dataSeries[,2])),2]<-1
dataSeries[which(is.na(dataSeries[,2])),2]<-0
#for (n in 2:(dim(dataSeries)[1])) if (dataSeries[n,2]>0) dataSeries[n,2]<-as.numeric(dataSeries[(n-1),2])+1
ncounter<-0
if (dataSeries[1,2]>0) ncounter<-1 # to avoid [0,]
#for (n in 2:(dim(dataSeries)[1])) if (dataSeries[n,2]>0) ncounter<-ncounter+1 else { if (dataSeries[n,2]>0) dataSeries[n,2]<-ncounter;ncounter<-0 }
for (n in 2:(dim(dataSeries)[1])) if (dataSeries[n,2]>0) ncounter<-ncounter+1 else { if (dataSeries[n-1,2]>0) dataSeries[n-1,2]<-ncounter;ncounter<-0 }
if (dataSeries[n,2]>0)  dataSeries[n,2]<-ncounter
if (dataSeries[n,2]==1)  dataSeries[n,2]<--1
for (m in 1:(dim(dataSeries)[1]-1)) if (dataSeries[m,2]==1 & dataSeries[m+1,2]==0) dataSeries[m,2]<--1
EndOfSegmentIndex<-which(dataSeries[,2]>1 | dataSeries[,2]<0)#get the indices where the segment sizes are stored
SegmentData<-dataSeries[EndOfSegmentIndex,]
SegmentData<-cbind.data.frame(StartDate=as.Date( as.numeric(as.Date(SegmentData[,1]))+1-abs(as.numeric(SegmentData[,2])) ,'1970-1-1'),EndDate=as.Date( SegmentData[,1]),ContinuousDays=abs(as.numeric(SegmentData[,2])))
SegmentData
}

Subgen<-function(select_flag)
{
#select_flag=c(0,1,1,0)
elements = c(',','i%d,')
selected_elements = elements[select_flag + 1]
format_str = paste(selected_elements,sep='',collapse='')
x<-list(substr(format_str,1,nchar(format_str)-1))
for (y in which(select_flag==1)) x<-c(x,y)
do.call(sprintf, x)
}

ndims<-function(x) length(dim(x))

BSXFUN<-function(op,x,y)
{# BSXFUN - apply a function to each element of 2 arrays with singleton expansion enabled
# original MATLAB code by Douglas M. Schwarz 2006
# translated to R by Jose Gama 2013
nd = max(ndims(x),ndims(y))
sx = dim(x)
if (length(sx)<nd) {
sx<-c(sx,rep(1,(nd-length(sx))))
x<-array(x,sx)
}
sy = dim(y)
if (length(sy)<nd) {
sy<-c(sy,rep(1,(nd-length(sy))))
y<-array(y,sy)
}
dz = sx != sy
dims = which(dz==TRUE)
num_dims = length(dims)

# Eliminate some simple cases.
if ((num_dims == 0) || (length(x) == 1) || (length(y) == 1)){
	z = op(x,y)
	return
}

# Check for dimensional compatibility of inputs, compute size and class of
# output array and allocate it.
if (!(all(sx[dz] == 1 || sy[dz] == 1))) stop('genop:argSizeError','Argument dimensions are not compatible.')

sz = apply(rbind(sx,sy),2,max)

#fix empty dimensions'
sz2 = apply(rbind(sx,sy),2,min)
if (any(sz2==0)) {
    sz[sz2==0]=0
    return(array(0,sz))
}

z1 = do.call(op,list(x[1],y[1]))
#if (is.logical(z1)) z =  matrix(FALSE,sz[1],sz[2]) else z = matrix(0,sz[1],sz[2])
if (is.logical(z1)) z =  array(FALSE,sz) else z = array(0,sz)

# Compute code strings representing the subscripts of x, y and z.
xsub = Subgen(sy != sz)
ysub = Subgen(sx != sz)
zsub = Subgen(dz)

# Generate the code.
indent = 2; # spaces per indent level
code_cells = vector('character',2*num_dims + 1)
for (i in 1:num_dims){
	code_cells[i] = sprintf('%*sfor (i%d in 1:sz[%d]){\n',indent*(i-1),'', dims[i], dims[i])
	code_cells[length(code_cells)-i+1] = sprintf('%*s}\n',indent*(i-1),'')
}
code_cells[num_dims+1] = sprintf('%*sz[%s] = do.call(op,list(x[%s],y[%s]))\n',indent*num_dims,'',zsub,xsub,ysub);
mycode = paste(code_cells, sep='',collapse='')
#mycode = paste('op<-"',op,'"\n',mycode, sep='',collapse='')
# Evaluate the code.
eval(parse( text=mycode ))
z
}

localwindow <- function(X, Y, DX, n){
# Index relevant to Local Window
Idx <- which(((X[n] - DX) <= X) & (X <= (X[n] + DX)))
# Calculate Local Nominal Data Reference Value
Y0  = median(Y[Idx])
# Calculate Local Scale of Natural Variation
S0  = 1.4826 * median(abs(Y[Idx] - Y0))
cbind(Y0=Y0, S0=S0)
}
     
smgauss <- function(X, V, DX){
# Prepare Xj and Xk
Xj  = matrix(X, length(X), length(X),byrow=T)
Xk  = matrix(X, length(X),length(X),byrow=F)
# Calculate Gaussian weight
Wjk = exp(-((Xj - Xk)/(2*DX))^2)
# Calculate Gaussian Filter
sWjk=apply(Wjk,1,sum)
m2=Wjk %*% V
G=m2 / sWjk
G
}

FindOutliersHampel <- function(X, Y, DX=NA, Th=NA, hampelAdaptive=FALSE, Threshold = NA)
{
# Hampel Filter.
# Chapters 1.4.2, 3.2.2 and 4.3.4 in Mining Imperfect Data: Dealing with 
# Contamination and Incomplete Records by Ronald K. Pearson.
# Ronald K. Pearson
# http://exploringdatablog.blogspot.com/2012/01/moving-window-filters-and-pracma.html
# Copyright (c) 2012:
# Michael Lindholm Nielsen
# http://www.mathworks.com/matlabcentral/fileexchange/34795-outlier-detection-and-removal--hampel-
N <- length(X)
SortX <- sort(X)
if (is.na(DX)) DX <- 3*median(SortX[2:N] - SortX[1:(N-1)])
if (is.na(Th)) Th <- 3
# Detect/Ignore NaN values in X and Y
IdxNaN  = which(is.na(X) | is.na(Y))
if (length(IdxNaN)>0){
X       = X[-IdxNaN]
Y       = Y[-IdxNaN]
}
# Calculation
# Preallocation
YY  = Y
I   = rep(FALSE,length(Y))
S0  = rep(NA,length(YY))
Y0  = S0
ADX = rep(DX, length(Y))

if (length(X) > 1){
if (hampelAdaptive==FALSE){ # standard
for (n in 1:length(Y)){
# Calculate Local Nominal Data Reference value
# and Local Scale of Natural Variation
tmp <- localwindow(X, Y, DX, n)
Y0[n] <- tmp[,'Y0']
S0[n] <- tmp[,'S0']
}
} else { # adaptive

# Preallocate
Y0Tmp   = S0
S0Tmp   = S0
DXTmp   = (1:length(S0)) * DX # Integer variation of Window Half Size

# Calculate Initial Guess of Optimal Parameters Y0, S0, ADX
for (n in 1:length(Y)){
# Setup/Reset temporary counter etc.
j   = 1
S0Rel   = Inf
while (S0Rel > Threshold){
# Calculate Local Nominal Data Reference value
# and Local Scale of Natural Variation using DXTmp window
tmp <- localwindow(X, Y, DXTmp[j], n)
Y0Tmp[j] <- tmp[,'Y0']
S0Tmp[j] <- tmp[,'S0']
# Calculate percent difference relative to previous value
if (j > 1) S0Rel <- abs((S0Tmp[j-1] - S0Tmp[j])/(S0Tmp[j-1] + S0Tmp[j])/2)
# Iterate counter
j   = j + 1
} # while
Y0[n]   = Y0Tmp[j - 2] # Local Nominal Data Reference value
S0[n]   = S0Tmp[j - 2] # Local Scale of Natural Variation
ADX[n] = DXTmp[j - 2]/DX  # Local Adapted Window size relative to DX
} # for
# Gaussian smoothing of relevant parameters
N <- length(SortX)
DX  = 2*median(SortX[2:N] - SortX[1:(N-1)])
ADX = smgauss(X, ADX, DX)
S0  = smgauss(X, S0, DX)
Y0  = smgauss(X, Y0, DX)
} # else
} # if

# Prepare Output
UB      = Y0 + Th*S0            # Save information about local scale
LB      = Y0 - Th*S0            # Save information about local scale
Idx     = which(abs(Y - Y0) > Th*S0)   # Index of possible outlier
YY[Idx] = Y0[Idx]              # Replace outliers with local median value
I[Idx]  = TRUE                 # Set Outlier detection
NO      = sum(I)               # Output number of detected outliers

# Reinsert NaN values detected at error checking stage
#if any(IdxNaN)
#    [YY, I, Y0, LB, UB, ADX]    = rescale(IdxNaN, YY, I, Y0, LB, UB, ADX);
cbind(YY=YY, I=I, Y0=Y0, LB=LB, UB=UB, ADX=ADX)
}

ESDstatistic <- function(dataSeries)
{# function to compute the test statistic ESD
# original code from: 
# NIST/SEMATECH e-Handbook of Statistical Methods, 2013
# http://www.itl.nist.gov/div898/handbook/
ares = abs(dataSeries - mean(dataSeries))/sd(dataSeries)
df = data.frame(y=dataSeries, ares)
r = max(df$ares)
list(r, df)
}

FindOutliersESDtest <- function(dataSeries,k=10,alpha=0.05)
{# generalized (extreme Studentized deviate) ESD test (Rosner 1983)
# original code from:
# NIST/SEMATECH e-Handbook of Statistical Methods, 2013
# http://www.itl.nist.gov/div898/handbook/
## Define values and vectors.
n = length(dataSeries)
alpha = 0.05
lam = c(1:k) # number of outliers, k
R = c(1:k)
## Compute test statistic until r=k values have been
## removed from the sample.
for (j in 1:k){
if(j==1){
rt = ESDstatistic(dataSeries)
R[j] = unlist(rt[1])
df = data.frame(rt[2])
newdf = df[df$ares!=max(df$ares),]}
else if(j!=1){
rt = ESDstatistic(newdf$y)
R[j] = unlist(rt[1])
df = data.frame(rt[2])
newdf = df[df$ares!=max(df$ares),]}
## Compute critical value.
p = 1 - alpha/(2*(n-j+1))
t = qt(p,(n-j-1))
lam[j] = t*(n-j) / sqrt((n-j-1+t**2)*(n-j+1))
}
## return results.
return( data.frame(No.Outliers=c(1:k),Test.Stat=R,Critical.Val=lam))
}

TietjenMoore <- function(dataSeries,k)
{# compute statistic Tietjen-Moore
# original code from: 
# NIST/SEMATECH e-Handbook of Statistical Methods, 2013
# http://www.itl.nist.gov/div898/handbook/
n = length(dataSeries)
## Compute the absolute residuals.
r = abs(dataSeries - mean(dataSeries))
## Sort data according to size of residual.
df = data.frame(dataSeries,r)
dfs = df[order(df$r),]
## Create a subset of the data without the largest k values.
klarge = c((n-k+1):n)
subdataSeries = dfs$dataSeries[-klarge]
## Compute the sums of squares.
ksub = (subdataSeries - mean(subdataSeries))**2
all = (df$dataSeries - mean(df$dataSeries))**2
## Compute the test statistic.
sum(ksub)/sum(all)
}

FindOutliersTietjenMooreTest <- function(dataSeries,k,alpha=0.05){
# Tietjen-Moore test
# original code from: 
# NIST/SEMATECH e-Handbook of Statistical Methods, 2013
# http://www.itl.nist.gov/div898/handbook/
ek <- TietjenMoore(dataSeries,k)
## Compute critical value based on simulation.
test = c(1:10000)
for (i in 1:10000){
dataSeriesdataSeries = rnorm(length(dataSeries))
test[i] = TietjenMoore(dataSeriesdataSeries,k)}
Talpha=quantile(test,alpha)
list(T=ek,Talpha=Talpha)
}

FindOutliersGrubbsTwosided <- function(dataSeries, alpha=0.05, iterative=TRUE)
{# returns the position of the estimated outliers using Grubb's iterative test for outliers
N <- length(dataSeries)
dataSeriesNA <- dataSeries
G <- max(abs(dataSeries - mean(dataSeries)))/sd(dataSeries)
posG <- which((abs(dataSeries - mean(dataSeries)))/sd(dataSeries)==G)[1] # locate the position of G
tcv=qt((1-alpha)/(2*N),N-2) # t statistic critical value
Gcv=((N-1)/sqrt(N))*sqrt(tcv^2/(N-2+tcv^2)) # G statistic critical value
if (G>Gcv) retG <- G else retG <- c()
if (!iterative) return(cbind(G=G,p=Gcv,posG=posG))
dataSeriesNA[posG] <- NA
dataSeries <- dataSeries[-posG] # remove the outlier G
retP <- Gcv
while((length(G)>0) & (N>=6))
{
N <- length(dataSeries)
G <- max(abs(dataSeries - mean(dataSeries)))/sd(dataSeries)
tcv=qt((1-alpha)/(2*N),N-2) # t statistic critical value
Gcv=((N-1)/sqrt(N))*sqrt(tcv^2/(N-2+tcv^2)) # G statistic critical value
if (length(G)>0) if (G>Gcv) {
retG <- c(retG,G)
retP <- c(retP,Gcv)
}
posG2 <- which((abs(dataSeries - mean(dataSeries)))/sd(dataSeries)==G)[1] # locate the position of G
posG3 <- which(dataSeries[posG2] == dataSeriesNA)[1] # locate the position of G in the original series
dataSeries <- dataSeries[-posG2] # remove the outlier G
dataSeriesNA[posG3] <- NA
if (length(G)>0) if (G>Gcv) posG <- c(posG,posG3)
}
cbind(G=retG,p=retP,posG=posG)
}

FindOutliersGrubbsOnesidedMax <- function(dataSeries, alpha=0.05, iterative=TRUE)
{# returns the position of the estimated outliers using Grubb's iterative test for outliers
N <- length(dataSeries)
dataSeriesNA <- dataSeries
G <- (max(dataSeries) - mean(dataSeries))/sd(dataSeries)
posG <- which((abs(dataSeries - mean(dataSeries)))/sd(dataSeries)==G)[1] # locate the position of G
tcv=qt((1-alpha)/(N),N-2) # t statistic critical value
Gcv=((N-1)/sqrt(N))*sqrt(tcv^2/(N-2+tcv^2)) # G statistic critical value
if (G>Gcv) retG <- G else retG <- c()
if (!iterative) return(cbind(G=G,p=Gcv,posG=posG))
dataSeriesNA[posG] <- NA
dataSeries <- dataSeries[-posG] # remove the outlier G
retP <- Gcv
while((length(G)>0) & (N>=6))
{
N <- length(dataSeries)
G <- (max(dataSeries) - mean(dataSeries))/sd(dataSeries)
tcv=qt((1-alpha)/(N),N-2) # t statistic critical value
Gcv=((N-1)/sqrt(N))*sqrt(tcv^2/(N-2+tcv^2)) # G statistic critical value
if (length(G)>0) if (G>Gcv) {
retG <- c(retG,G)
retP <- c(retP,Gcv)
}
posG2 <- which((abs(dataSeries - mean(dataSeries)))/sd(dataSeries)==G)[1] # locate the position of G
posG3 <- which(dataSeries[posG2] == dataSeriesNA)[1] # locate the position of G in the original series
dataSeries <- dataSeries[-posG2] # remove the outlier G
dataSeriesNA[posG3] <- NA
if (length(G)>0) if (G>Gcv) posG <- c(posG,posG3)
}
cbind(G=retG,p=retP,posG=posG)
}

FindOutliersGrubbsOnesidedMin <- function(dataSeries, alpha=0.05, iterative=TRUE)
{# returns the position of the estimated outliers using Grubb's iterative test for outliers
N <- length(dataSeries)
dataSeriesNA <- dataSeries
G <- (mean(dataSeries) - min(dataSeries))/sd(dataSeries)
posG <- which((abs(dataSeries - mean(dataSeries)))/sd(dataSeries)==G)[1] # locate the position of G
tcv=qt((1-alpha)/(N),N-2) # t statistic critical value
Gcv=((N-1)/sqrt(N))*sqrt(tcv^2/(N-2+tcv^2)) # G statistic critical value
if (G>Gcv) retG <- G else retG <- c()
if (!iterative) return(cbind(G=G,p=Gcv,posG=posG))
dataSeriesNA[posG] <- NA
dataSeries <- dataSeries[-posG] # remove the outlier G
retP <- Gcv
while((length(G)>0) & (N>=6))
{
N <- length(dataSeries)
G <- (mean(dataSeries) - min(dataSeries))/sd(dataSeries)
tcv=qt((1-alpha)/(N),N-2) # t statistic critical value
Gcv=((N-1)/sqrt(N))*sqrt(tcv^2/(N-2+tcv^2)) # G statistic critical value
if (length(G)>0) if (G>Gcv) {
retG <- c(retG,G)
retP <- c(retP,Gcv)
}
posG2 <- which((abs(dataSeries - mean(dataSeries)))/sd(dataSeries)==G)[1] # locate the position of G
posG3 <- which(dataSeries[posG2] == dataSeriesNA)[1] # locate the position of G in the original series
dataSeries <- dataSeries[-posG2] # remove the outlier G
dataSeriesNA[posG3] <- NA
if (length(G)>0) if (G>Gcv) posG <- c(posG,posG3)
}
cbind(G=retG,p=retP,posG=posG)
}

FindOutliersZscore <- function(dataSeries, coef=2.5)
{# returns the position of the values outside the allowed range by a criteria based on the z score, abs(Z)>coef
zscore <- scale(dataSeries)
which((abs(zscore)>coef))
}

FindOutliersModifiedZscore <- function(dataSeries, coef=3.5)
{# returns the position of the values outside the allowed range by a criteria based on the modified z score, abs(Z)>coef
zscore <- 0.6745*(dataSeries-median(dataSeries))/mad(dataSeries)
which((abs(zscore)>coef))
}

TrimmedMean <- function(dataSeries, percentTrim=0.1)
{ # trimmed mean
N <- length(dataSeries)
NP <- round(N * percentTrim,1)
R<- N-2*NP
dataSeries <- dataSeries[-c(1:floor(NP),(N-floor(NP)+1):N)]
if (NP-floor(NP)!=0) dataSeries[c(1:floor(NP),(length(dataSeries)-floor(NP)+1))] <- dataSeries[c(1:floor(NP),(length(dataSeries)-floor(NP)+1))]*(NP-floor(NP))
T <- sum(dataSeries)/(R)
T
}

FindOutliersQuant <- function(dataSeries, coef=1.5)
{# returns the position of the values outside the allowed range by a criteria based on quantiles, q25-coef*(q75-q25)<x<q75+coef*(q75-q25)
q25 <- quantile(dataSeries,.25, na.rm =T)
q75 <- quantile(dataSeries,.75, na.rm =T)
which((dataSeries<(q25-coef*(q75-q25)))|(dataSeries>(q75+coef*(q75-q25))))
}

FindOutliersSD <- function(dataSeries, coef=3)
{# returns the position of the values outside the allowed range by a criteria based on the standard deviation, mean-coef*sd<x<mean+coef*sd
seriesSD <- sd(dataSeries)
seriesMean <- mean(dataSeries)
which((dataSeries<(seriesMean-coef*seriesSD))|(dataSeries>(seriesMean+coef*seriesSD)))
}

FindOutliersTrimmedMeans <- function(dataSeries, percentTrim=0.1,coef=3)
{# returns the positions of the values outside the allowed range by a criteria based on trimmed means
seriesTM <- TrimmedMean(dataSeries, percentTrim)
seriesSD <- sd(dataSeries)
which((dataSeries<(seriesTM-coef*seriesSD))|(dataSeries>(seriesTM+coef*seriesSD)))
}

FindOutliersMAD <- function(dataSeries, coef=3)
{# returns the position of the values outside the allowed range by a criteria based on the absolute deviation around the median (MAD), median-coef*MAD<x<median+coef*MAD
seriesMAD <- mad(dataSeries)
seriesMedian <- median(dataSeries)
which((dataSeries<(seriesMedian-coef*seriesMAD))|(dataSeries>(seriesMedian+coef*seriesMAD)))
}

CompleteMissingValuesDailyMean <- function(dataSeries, windowBefore=1,windowAfter=1, missingValue=-9999)
{ # complete missing values on a time series of daily values by using the average of 1 or more days before or after
# the missing values remains unchanged if there is one or more missing values in the chosen window of time for calculating the mean
N <- dim(dataSeries)[1]
dataSeries <- FillDailyGapsWithSomeValue(dataSeries,dataSeries[1,1], dataSeries[N,1], missingValue=missingValue)
dataSeries[which(dataSeries[,2]==missingValue),2] <- NA # force it to be NA or it will cause an error on sapply if missingValue=NA
rangeData <- (windowBefore+1):(N-windowAfter)
dataSeries[,2] <- sapply(rangeData, function(x) {
currentX <-ifelse(windowBefore==0,1,windowBefore+1)
ifelse(is.na(dataSeries[x,2]),mean(dataSeries[c((x-windowBefore):(x+windowAfter))[-currentX],2]),dataSeries[x,2])})
dataSeries
}

RankDifferenceTest <- function(dataSeries)
{# Rank Difference (non-parametric test for randomness)
# Trend 1.0.2 User Guide, chapter 4.2.11 Rank Difference Test, pp. 21
# http://www.toolkit.net.au/Tools/TREND/documentation
N <- length(dataSeries)
ranksX <- rank(dataSeries)
ranksDiff <- abs(diff(ranksX))
U <- sum(ranksDiff)
m <- (N+1)*(N-1)/3
sigma <- (N-2)*(N+1)*(4*N-7)/90
z <- abs(U-m)/sqrt(sigma)
list(m=m,sigma=sigma,z=z)
}

TurningPointsTest <- function(dataSeries)
{# Turning Points (non-parametric test for randomness)
# Trend 1.0.2 User Guide, chapter 4.2.10 Turning Points Test, pp. 21
N <- length(dataSeries)
x2 <- sapply(2:(N-1),function(n) ifelse(((dataSeries[n]>dataSeries[n-1]) & (dataSeries[n]>dataSeries[n+1])) | ((dataSeries[n]<dataSeries[n-1]) & (dataSeries[n]<dataSeries[n+1])),1,0))
mS <- sum(x2)
m <- 2*(N-2)/3
sigma <- (16*N-29)/90
z <- abs(mS-m)/sqrt(sigma)
list(m=m,sigma=sigma,z=z)
}

MedianCrossingTest <- function(dataSeries)
{# Median Crossing (non-parametric test for randomness)
# Trend 1.0.2 User Guide, chapter 4.2.9 Median Crossing Test, pp. 21
medianX <- median(dataSeries)
x2 <- ifelse(dataSeries>medianX,1,0)
N <- length(dataSeries)
z <- sapply(2:N, function(n) ifelse(x2[n]==x2[n-1],0,1))
z
}

RankSumTest <- function(dataSeries, period1)
{# Rank-Sum (non-parametric test for difference in median from two data periods)
# Trend 1.0.2 User Guide, chapter 4.2.7 Rank-Sum Test, pp. 20
N <- length(dataSeries)
period2 <- N - period1
ranksX <- rank(dataSeries)
S <- sum(ranksX[(period2+1):N])
m <- period1 * (N+1)/2
sigma <- sqrt(period1*period2*(N+1)/12)
z <- ifelse((S>m),(S-0.5-m)/sigma,(S+0.5-m)/sigma)
list(S=S,m=m,sigma=sigma,z=z)
}

CreateReferenceSeriesFromFilesMeanCorrelations<-function(vFiles,commonPeriod=NA,refSeriesFile=NA,wholePeriod=FALSE, useDiff=TRUE,retInfo=FALSE)
{ # creates a reference series from two or more series
# same as AnClim 
tmp <- FindCommonPeriod(vFiles,TRUE)
commonPeriodCalc <- tmp$commonPeriod
allMin <- tmp$allMin
allMax <- tmp$allMax
if (any(is.na(commonPeriod))) commonPeriod <- tmp$commonPeriod else {
if (commonPeriodCalc[1]>commonPeriod[1]) commonPeriod[1] <- commonPeriodCalc[1]
if (commonPeriodCalc[2]<commonPeriod[2]) commonPeriod[2] <- commonPeriodCalc[2]
}
if (wholePeriod) commonPeriod <- c(allMin,allMax)
CommonPeriodSeq <- commonPeriod[1]:commonPeriod[2]
N <- length(vFiles)
fmN<-read.delim(vFiles[N],sep='',header=FALSE)
fmN<-fmN[which(fmN[,1] %in% CommonPeriodSeq),]
meanfmN <- apply(fmN,2,mean)[-1]
allMeans <- matrix(0,dim(fmN)[2]-1,N)
allMeans[,N] <- meanfmN
allCorr <- matrix(0,dim(fmN)[2]-1,N-1)
allWeights <- matrix(0,dim(fmN)[2]-1,N-1)
w2=0
for (n in 1:(N-1)){
fm1<-read.delim(vFiles[n],sep='',header=FALSE)
fm1 <- fm1[which(fm1[,1] %in% CommonPeriodSeq),]
meanfm1 <- apply(fm1,2,mean)[-1]
allMeans[,n] <- meanfm1
allCorr[,n] <- sapply(2:(dim(fm1)[2]), function(x) cor(fm1[,x],fmN[,x])) 
allWeights[,n] <- sapply(2:(dim(fm1)[2]), function(x) sum( cor(fm1[,x],fmN[,x])^2 ))
mfm1=matrix(c(0,meanfm1),dim(fm1)[1],dim(fm1)[2],byrow=T)
mfmN=matrix(c(0,meanfmN),dim(fm1)[1],dim(fm1)[2],byrow=T)
if (useDiff) mfm3 <- fm1-mfm1+mfmN else mfm3 <- fm1*mfm1/mfmN
mw=matrix(c(1,allCorr[,n]),dim(fm1)[1],dim(fm1)[2],byrow=T)^2
w2=w2+ mfm3 * mw
}
wsum <- apply(allCorr,1, function(x) sum(x^2))
w2=w2/matrix(c(1,wsum),dim(fm1)[1],dim(fm1)[2],byrow=T)
w2 <- round(w2,1)
 cat('Weighted-average calculated differences\nSources for the reference series:\n\n',vFiles)
 for (n in 1:(dim(fm1)[2]-1)){
 cat('Series',n,'files 1 -',N,', years',commonPeriod[1],'-',commonPeriod[2],'\n')
 cat('Mean:',allMeans[n,],'\n')
 cat('Weights/Correlations:',allWeights[n,], allCorr[n,],'\n')
 }
 if (!is.na(refSeriesFile)) write.table(fm1,refSeriesFile,w2, row.names = FALSE, col.names = FALSE)
 if (!retInfo) return(w2) else return(list(refSeries=w2, meansSeries=allMeans, corrSeries=allCorr, weightsSeries=allWeights))
}

CreateReferenceSeriesFromFilesMeanCorrelationsTwoseries<-function(vFiles,commonPeriod=NA,refSeriesFile=NA,wholePeriod=FALSE, useDiff=TRUE,retInfo=FALSE)
{ # creates a reference series from two series
# same as AnClim "Create reference series"
tmp <- FindCommonPeriod(vFiles,TRUE)
commonPeriodCalc <- tmp$commonPeriod
allMin <- tmp$allMin
allMax <- tmp$allMax
if (any(is.na(commonPeriod))) commonPeriod <- tmp$commonPeriod else {
if (commonPeriodCalc[1]>commonPeriod[1]) commonPeriod[1] <- commonPeriodCalc[1]
if (commonPeriodCalc[2]<commonPeriod[2]) commonPeriod[2] <- commonPeriodCalc[2]
}
if (wholePeriod) commonPeriod <- c(allMin,allMax)
CommonPeriodSeq <- commonPeriod[1]:commonPeriod[2]
fm1<-read.delim(vFiles[1],sep='',header=FALSE)
fm1<-fm1[which(fm1[,1] %in% CommonPeriodSeq),]
meanfm1 <- apply(fm1,2,mean)[-1]
w <- meanfm1
N <- length(vFiles)
allMeans <- matrix(0,dim(fm1)[2]-1,N)
allMeans[,1] <- meanfm1
allCorr <- matrix(0,dim(fm1)[2]-1,N-1)
allWeights <- matrix(0,dim(fm1)[2]-1,N-1)
for (n in 2:N){
fm2<-read.delim(vFiles[n],sep='',header=FALSE)
fm2<-fm2[which(fm2[,1] %in% CommonPeriodSeq),]
meanfm2 <- apply(fm2,2,mean)[-1]
allMeans[,n] <- meanfm2
w <- w-meanfm2
allCorr[,n-1] <- sapply(2:dim(fm1)[2], function(x) cor(fm1[,x],fm2[,x])) 
allWeights[,n-1] <- sapply(2:dim(fm1)[2], function(x) sum( cor(fm1[,x],fm2[,x])^2 ))
}
w <- round(w,1)
if (useDiff) refScalc <- fm1-matrix(c(0,w),dim(fm1)[1],dim(fm1)[2],byrow=T) else refScalc <- fm1/matrix(c(0,w),dim(fm1)[1],dim(fm1)[2],byrow=T)
cat('Weighted-average calculated differences\nSources for the reference series:\n\n',vFiles)
for (n in 1:(dim(fm1)[2]-1)){
cat('Series',n,'files 1 -',N,', years',commonPeriod[1],'-',commonPeriod[2],'\n')
cat('Mean:',allMeans[n,],'\n')
cat('Weights/Correlations:',allWeights[n,], allCorr[n,],'\n')
}
if (!is.na(refSeriesFile)) write.table(fm1,refSeriesFile,refScalc, row.names = FALSE, col.names = FALSE)
if (!retInfo) return(refScalc) else return(list(refSeries=refScalc, meansSeries=allMeans, corrSeries=allCorr, weightsSeries=allWeights))
}

CumulativeDeviations <- function(dataSeries) abs(cumsum(dataSeries-mean(dataSeries))) / sqrt(sum((dataSeries-mean(dataSeries))^2)/length(dataSeries))# returns the cumulative deviations test (parametric test for step jump in mean) value Q
CumulativeDeviationsQR <- function(dataSeries)
{ # returns the cumulative deviations test (parametric test for step jump in mean) value Q, range R and Q/sqrt(n)
# Trend 1.0.2 User Guide, chapter 4.2.4, 4.2.5 Cumulative Deviation Test, pp. 18
Q <- (cumsum(dataSeries-mean(dataSeries))) / sqrt(sum((dataSeries-mean(dataSeries))^2)/length(dataSeries))
list(Q=abs(Q),R=abs(Q-min(Q)),QsqrtN=abs(Q)/sqrt(length(dataSeries)))
}

DistributionFreeCUSUM <- function(dataSeries)
{ # returns the Distribution Free CUSUM (non-parametric test for step jump in mean)
# Trend 1.0.2 User Guide, chapter 4.2.4, Distribution Free CUSUM Test, pp. 18
s <- sign(dataSeries-median(dataSeries))
vK <- cumsum(s)
absVK <- abs(vK)
Q <- max(absVK)
Q
}

# cumulativeDeviationsTestQR <- function(x)
# { # returns the cumulative deviations test value Q, range R and Q/sqrt(n)
# n=length(x)
# xMean <- x-mean(x)
# xMean2 <- xMean^2
# sK <- cumsum(xMean)
# Dxr <- sqrt(sum(xMean2)/n)
# sK2 <- sK / Dxr
# absSK2 <- abs(sK2)
# Q <- max(absSK2)
# list(Q=Q,R=Q-min(absSK2),QsqrtN=Q/sqrt(n))
# }

FillDailyGapsWithSomeValue<-function(dataYearSeries,FromDate,ToDate, missingValue=-9999)
{#fill missing days with date+missingValue when the timse series is literally missing these days
if (is.character(FromDate)) FromDate<-as.Date(FromDate, format="%Y-%m-%d")
if (is.character(ToDate)) ToDate<-as.Date(ToDate, format="%Y-%m-%d")
dataYearSeries<-ValuesBetween2years(dataYearSeries,FromDate,ToDate)
AllDaysMeasuredBetween1stLastDays<-dataYearSeries[,1]
firstDay<-min(AllDaysMeasuredBetween1stLastDays, na.rm =T)
lastDay<-max(AllDaysMeasuredBetween1stLastDays, na.rm =T)
AllDaysBetween1stLastDays<-seq(firstDay,lastDay,by="1 day")
tseries2<-data.frame(cbind(as.character(AllDaysBetween1stLastDays),NA))
#too slow
#for (n in 1:dim(tseries2)[1]) { z<-which(dataYearSeries[,1]==as.Date(tseries2[n,1])); if (length(z)>0) tseries2[n,2]<-dataYearSeries[z,1] }
MissingDates<-as.character(tseries2[,1]) %w/o% as.character(dataYearSeries[,1])
if (length(MissingDates)==0) return (dataYearSeries)
MissingDates<-data.frame(cbind(MissingDates,NA))
colnames(MissingDates)<-colnames(dataYearSeries);
#colnames(dataYearSeries)<-colnames(MissingDates);
##rownames(dataYearSeries)<-rownames(MissingDates);
dataYearSeries<-rbind.data.frame(dataYearSeries,as.data.frame(MissingDates))
if (!is.na(missingValue)) dataYearSeries[which(is.na(dataYearSeries[,2])),2] <- missingValue # replace NA with missingValue
dataYearSeries[,2]<-as.numeric(dataYearSeries[,2])
dataYearSeries[order(dataYearSeries[,1]),]
}

FillYearlyGapsWithSomeValue<-function(dataYearSeries,FromYear=min(dataYearSeries[,1]),ToYear=max(dataYearSeries[,1]), missingValue=-9999)
{#fill missing years with year+missingValue when the timse series is literally missing these years
  MissingDates <- FromYear:ToYear  %w/o% dataYearSeries[,1]
  MissingDates<-data.frame(cbind(MissingDates,NA))
  colnames(MissingDates) <- colnames(dataYearSeries)
  dataYearSeries<-rbind.data.frame(dataYearSeries,as.data.frame(MissingDates))
  if (!is.na(missingValue)) dataYearSeries[which(is.na(dataYearSeries[,2])),2] <- missingValue # replace NA with missingValue
  dataYearSeries[order(dataYearSeries[,1]),]
}

ValuesBetween2Dates<-function(dataYearSeries,FromDate,ToDate)
{# return all the values between two dates, inclusively
dtfrm2<-dataYearSeries[which(dataYearSeries[,1]<=ToDate),]
dtfrm2<-dtfrm2[which(dtfrm2[,1]>=FromDate),]
dtfrm2
}

ValuesBetween2years<-function(dataYearSeries,intFromYear,intToYear)
{# return all the values between 2 years, that is, the 1st and the last days of a year range
if (is.numeric(dataYearSeries[1,1])) 
{
tmp<-as.character(dataYearSeries[,1])
tmp<-as.Date(paste(tmp, '-01-01',sep=''), format="%Y-%m-%d")
tmp<-cbind.data.frame(tmp, dataYearSeries[,2:dim(dataYearSeries)[2]])
} else tmp<-dataYearSeries
dtfrm2<-tmp[which(tmp[,1]<paste(intToYear+1,'01-01',sep='-')),]
dtfrm2<-dtfrm2[which(dtfrm2[,1]>=paste(intFromYear,'01-01',sep='-')),]
if (is.numeric(dataYearSeries[1,1])) dtfrm2[,1]<-as.numeric(format(as.Date(dtfrm2[,1]), "%Y") )
dtfrm2
}

FindCommonPeriod<-function(vFiles, returnMaxMin=FALSE)
{ # returns the common period for multiple files having the year in the first column
CommonPeriod <- c(Inf,-Inf)
vMin <- Inf
vMax <- -Inf
for (n in 1:length(vFiles)){
fm1<-read.delim(vFiles[n],sep='',header=FALSE)
CommonPeriod <- c(max(min(fm1[,1]),min(CommonPeriod)),min(max(fm1[,1]),max(CommonPeriod)))
if (min(fm1[,1])<vMin) vMin <- min(fm1[,1])
if (max(fm1[,1])>vMax) vMax <- max(fm1[,1])
}
if (returnMaxMin) return (list(commonPeriod=CommonPeriod,allMin=vMin,allMax=vMax)) else return(CommonPeriod)
}

# FillYearlyGapsWithNA<-function(dataYearSeries,fromYear,toYear)
# { # fill missing years with year + 12 months of NA
# y <- fromYear:toYear
# tmp <- data.frame(years=y,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
# colnames(tmp) <- c('years',paste('m',1:12,sep=''))
# tmp[which(tmp[,1] %in% dataYearSeries[,1]),] <- dataYearSeries
# tmp
# }

CreateReferenceSeriesFromFilesMean<-function(vFiles,commonPeriod=NA,refSeriesFile=NA,wholePeriod=FALSE,deviationsFlag=FALSE)
{ # creates a reference series from two or more series
tmp <- FindCommonPeriod(vFiles,TRUE)
commonPeriodCalc <- tmp$commonPeriod
allMin <- tmp$allMin
allMax <- tmp$allMax
if (any(is.na(commonPeriod))) commonPeriod <- tmp$commonPeriod else {
if (commonPeriodCalc[1]>commonPeriod[1]) commonPeriod[1] <- commonPeriodCalc[1]
if (commonPeriodCalc[2]<commonPeriod[2]) commonPeriod[2] <- commonPeriodCalc[2]
}
if (wholePeriod) commonPeriod <- c(allMin,allMax)
CommonPeriodSeq <- commonPeriod[1]:commonPeriod[2]
fm1<-read.delim(vFiles[1],sep='',header=FALSE)
if (deviationsFlag){
meanfm1 <- apply(fm1,2,mean)[-1]
fm1[,-1]<-fm1[,-1]-matrix(meanfm1,dim(fm1)[1],ncol=12,byrow=T)
}
if (wholePeriod) fm1<- FillYearlyGapsWithSomeValue(fm1,allMin,allMax) else fm1<-fm1[which(fm1[,1] %in% CommonPeriodSeq),]
N <- length(vFiles)
for (n in 2:N){
fm2<-read.delim(vFiles[n],sep='',header=FALSE)
fm2<-fm2[which(fm2[,1] %in% CommonPeriodSeq),]
if (deviationsFlag){
meanfm2 <- apply(fm2,2,mean)[-1]
fm2[,-1]<-fm2[,-1] - matrix(meanfm2,dim(fm1)[1],ncol=12,byrow=T)
}
if (wholePeriod) fm2 <- FillYearlyGapsWithSomeValue(fm2,allMin,allMax)
fm1 <- fm1+fm2
}
fm1<-fm1/N
colnames(fm1) <- c('years',paste('m',1:12,sep=''))
if (wholePeriod) fm1[which(is.na(fm1)) %% dim(fm1)[1] ,-1] <- -999.0
if (is.na(refSeriesFile)) return(fm1) else write.table(fm1,refSeriesFile,row.names = FALSE, col.names = FALSE)
}

PlotSeriesAndRefSimple <- function(dataYearSeries=NA,refYearSeries=NA)
{ # plot a series and one or more reference series
# similar to Anclim's Analyzing Plot Series
mn1 <- min(dataYearSeries[,2],refYearSeries[,-1])
mx1 <- max(dataYearSeries[,2],refYearSeries[,-1])
yL <- c(mn1,mx1)
plot(dataYearSeries[,1],refYearSeries[,2],col=1,type='l',ylim=yL);par(new=TRUE)
for (n in 2:dim(refYearSeries)[2]){
plot(dataYearSeries[,1],dataYearSeries[,n],col=n,lty=n,type='l',ylim=yL)
if (n<dim(refYearSeries)[2]) par(new=TRUE)
}
}

PlotSeriesDifference <- function(dataYearSeries=NA,refYearSeries=NA, diffFlag=TRUE)
{ # plot the difference between 2 series (series - reference series)
# similar to Anclim's Plot of Differences
if (diffFlag) tTemp <- dataYearSeries[,2]-refYearSeries[,2] else tTemp <- dataYearSeries[,2] / refYearSeries[,2]
yL <- c(min(tTemp),max(tTemp))
plot(dataYearSeries[,1],tTemp,col='blue',type='l',ylim=yL)
}

SNHTrelative<-function(dataSeries=NA,refSeries=NA, diffFlag=TRUE)
{#calculates the SNHT relative, assumes no NAs
if (diffFlag) SNHTabsolute(dataSeries-refSeries) else SNHTabsolute(dataSeries/refSeries)
}

SNHTrelativeII<-function(dataSeries=NA,refSeries=NA, diffFlag=TRUE)
{#calculates the SNHT relative, standard deviation different from 1, assumes no NAs
if (diffFlag) SNHTabsoluteII(dataSeries-refSeries) else SNHTabsoluteII(dataSeries/refSeries)
}

SNHTrelative2SDs<-function(dataSeries=NA,refSeries=NA, diffFlag=TRUE)
{#calculates the SNHT relative, two different standard deviations, assumes no NAs
if (diffFlag) SNHTabsolute2SDs(dataSeries-refSeries) else SNHTabsolute2SDs(dataSeries/refSeries)
}

SNHTrelativeDoubleShift<-function(dataSeries=NA,refSeries=NA, diffFlag=TRUE)
{#calculates the SNHT relative, Double shift of the mean level, assumes no NAs
if (diffFlag) SNHTabsoluteDoubleShift(dataSeries-refSeries) else SNHTabsoluteDoubleShift(dataSeries/refSeries)
}

SNHTrelativeTrend<-function(dataSeries=NA,refSeries=NA, diffFlag=TRUE)
{#calculates the SNHT relative with trend in the mean level, assumes no NAs
if (diffFlag) SNHTabsoluteTrend(dataSeries-refSeries) else SNHTabsoluteTrend(dataSeries/refSeries)
}

Sqplot <- function(n){
# returns the number of elements for a plot to be nearly square
n1 <- round(sqrt(n))
if ((n1^2) >=n) return(c(n1,n1))
return(c(n1,n1+1))
}

PlotHomog<-function(dataYearSeries=NA,refYearSeries=NA,diffFlag=TRUE,homogenization=SNHTabsolute,levelSignificance=c(99,95),criticalValues=climtrends::SNHT.Critical.Values,posLegend="topright",rowbycol=NULL)
{# plot the graphics for many homogenization methods with an indication of the peak and significance levels
if (any(is.na(dataYearSeries))) stop('<<<dataYearSeries>>> must be a table or dataframe with years and values, no NAs')
if (any(is.na(refYearSeries))) refYearSeries <- NA
if (!is.na(refYearSeries)) if (diffFlag) dataYearSeries[,2] <- dataYearSeries[,2]-refYearSeries[,2] else dataYearSeries[,2] <- dataYearSeries[,2]/refYearSeries[,2]
if (is.function(homogenization)) mystrfun<-as.character(substitute(homogenization)) else stop('<<<homogenization>>> must be a valid function name')
if (mystrfun=='homogenization') {
    nn <- substitute(homogenization)
    i <- 1
    while(TRUE) {
      on <- do.call("substitute", list(as.name(nn), parent.frame(i)))
      if (on == nn) break;
      nn <- on
      i <- i + 1
    }
    mystrfun<-nn
  }
critSNHT<-as.numeric(rownames(criticalValues))
N <- dim(dataYearSeries)[1]
N2 <- dim(dataYearSeries)[2]
minD <- min(abs(critSNHT-N))
Nsnht <- as.character(critSNHT[which(minD==abs(critSNHT-N))])
if (length(Nsnht)>1) Nsnht <- Nsnht[2]
if (is.numeric(levelSignificance)) levelSignificance <- paste(levelSignificance,ifelse(mystrfun=='BuishandRangeTest','R%','%'),sep='')
c1 <- criticalValues[Nsnht, levelSignificance]
if (is.null(rowbycol)) rowbycol <- Sqplot(N2-1)
s1 <- matrix(NA,N,N2)
s1[,1] <- dataYearSeries[,1]
for (n in 1:(N2-1)){
stmp <- homogenization(dataYearSeries[,(1+n)])
print(stmp)
s2 <- NA
if (is.list(stmp)) {
s2 <- unlist(stmp[[2]])
stmp <- unlist(stmp[[1]])
}
s1[1:length(stmp),1+n] <- stmp
}
# plot the SNHT T
if (mystrfun %in% c('BuishandRangeTest','PettittTest')) peakPointFunction <- min else peakPointFunction <- max
par(mfrow=rowbycol,mar=c(1,1,1,1),oma=c(1, 3, 0,0))
for (n in 1:(N2-1)) {
plot(s1[,1], s1[,n+1] ,type='l')
abline(v=s1[,1],col='lightgrey',lty=3)
abline(v=(dataYearSeries[which(s1[,n+1]==peakPointFunction(s1[,n+1], na.rm = TRUE)),1]),col='red')
if (!is.na(s2)) abline(v=(dataYearSeries[s2[which(s1[,n+1]==peakPointFunction(s1[,n+1], na.rm = TRUE))],1]),col='red',lty=2)
if ((peakPointFunction(s1[,n+1], na.rm = TRUE))>=(max(c1)))
{
legend(posLegend,legend=levelSignificance, text.col = 'green',title='',lty=1:2,box.lty=0)
abline(h=c1[1],col='green')
abline(h=c1[2],col='green',lty=2)
} else if ((peakPointFunction(s1[,n+1], na.rm = TRUE))>=(min(c1)))
{
legend(posLegend,legend=min(levelSignificance), text.col = 'green',title='',lty=2,box.lty=0)
abline(h=min(c1),col='green',lty=2)
}
}
}

AdjHomogenization<-function(dataYearSeries=NA,yearShift=NA,nYears=20, diffFlag=TRUE)
{ # get the value for SNHT adjustment within a certain window
# diifFlag=TRUE means difference, FALSE means ratio
if (any(is.na(dataYearSeries))) stop('<<<dataYearSeries>>> must be a table or dataframe with years and values, no NAs')
if (is.na(yearShift)) stop('<<<yearShift>>> must be a scalar numeric')
if (length(yearShift)>1) {
stop('<<<yearShift>>> must be an integer scalar')}
if (all(!(yearShift %in% dataYearSeries[,1]))) stop('<<<yearShift>>> must be within the year range')
nShift<-which(dataYearSeries[,1]==yearShift) # get pos of year shift
# window
nMin <- nShift-nYears
if (nMin<1) nMin <- 1
nMax <- nShift+nYears-1
if (nMax>dim(dataYearSeries)[1]) nMax <- dim(dataYearSeries)[1]
#cat('*',nShift,nMax,nMin,(nShift-1),'*')
if (diffFlag) return(mean(dataYearSeries[(nShift):nMax,2]) - mean(dataYearSeries[nMin:(nShift-1),2])) else return(mean(dataYearSeries[(nShift):nMax,2]) / mean(dataYearSeries[nMin:(nShift-1),2]))
}

GetInfoHomogenization<-function(dataYearSeries=NA,refYearSeries=NA,nYears=20, diffFlag=TRUE,returnData=FALSE, homogenization=SNHTabsolute,levelSignificance=c(99,95),criticalValues=climtrends::SNHT.Critical.Values)
{ # calculate SNHT or other method of homogenization and year of change, To value, adjust value, etc
#if (is.null(criticalValues)) criticalValues <-  get("SNHT.Critical.Values", envir  = environment())

if (any(is.na(dataYearSeries))) stop('<<<dataYearSeries>>> must be a table or dataframe with years and values, no NAs')
if (any(is.na(refYearSeries))) refYearSeries <- NA
if (!is.na(refYearSeries)) if (diffFlag) dataYearSeries[,2] <- dataYearSeries[,2]-refYearSeries[,2] else dataYearSeries[,2] <- dataYearSeries[,2]/refYearSeries[,2]
if (is.function(homogenization)) mystrfun<-as.character(substitute(homogenization)) else stop('<<<homogenization>>> must be a valid function name')
if (mystrfun=='homogenization') {
    nn <- substitute(homogenization)
    i <- 1
    while(TRUE) {
      on <- do.call("substitute", list(as.name(nn), parent.frame(i)))
      if (on == nn) break;
      nn <- on
      i <- i + 1
    }
    mystrfun<-nn
  }
N <- dim(dataYearSeries)[1]
critSNHT<-as.numeric(rownames(criticalValues))
minD <- min(abs(critSNHT-N))
Nsnht <- as.character(critSNHT[which(minD==abs(critSNHT-N))])
if (length(Nsnht)>1) Nsnht <- Nsnht[2]
if (is.numeric(levelSignificance)) levelSignificance <- paste(levelSignificance,ifelse(mystrfun=='BuishandRangeTest','R%','%'),sep='')
#cat('\t*******\t',minD,levelSignificance,'\n')
c1 <- criticalValues[Nsnht, levelSignificance]
#c5 <- criticalValues[Nsnht,"95%"]
s1 <- homogenization(dataYearSeries[,2])
# cases when there is more than 1 output (2 shifts or trend)
s2 <- NULL
if (is.list(s1)) {s2 <- unlist(s1[[2]]); s1 <- unlist(s1[[1]])}
# max or min depends on the test
if ((mystrfun=='PettittTest')|(mystrfun=='BuishandRangeTest')) ToV <- min(s1) else ToV <- max(s1)
# in case of unknown tests also with negative peaks
if ((ToV==0) & (NegCurve(s1)==-1)) ToV <- min(s1)
yearChange <- dataYearSeries[which(s1 == ToV),1]+1
if (length(yearChange)>1){# if there is a tie between 2 years, use the Z statistic to choose one
Zstat <- scale(dataYearSeries[,2])*NegCurve(s1)
  z2 <- Zstat[which(s1 == ToV)]
  z3 <- order(z2)
  yearChange <- yearChange[z3]
#   print(s1)
#   print(z2)
#   print(z3)
#   print(ToV)
#   print(yearChange)
  yearChange <- yearChange[1]
  #plot(dataYearSeries[,1],s1,type='l')
  #;par(new=TRUE);plot(dataYearSeries[,1],Zstat,type='l',col='grey');
}
#cat(yearChange,ToV,'\n')
# the year might be not within the range
if (all(!(yearChange %in% dataYearSeries[,1]))) AdjustValue <- 0 else AdjustValue <- AdjHomogenization(dataYearSeries,yearChange,nYears,diffFlag)
SNHTret <- list(ToValue=ToV, yearChange=yearChange,criticalFound=c1,AdjustValue=AdjustValue)
if (returnData) SNHTret$SNHTdata <- s1
if (!is.null(s2)) SNHTret$endyear <- dataYearSeries[s2[(which(s1 == ToV))],1]+1
SNHTret
}

LevelSignifMarks <- function(dataSeries=NA, levelSignificance=c(99,95), userDefSymbol='<')
{# returns a certain number of marks depending on the level of significance achieved by the value statistic
if (is.na(dataSeries)) return('')
for (n in length(levelSignificance)) if (dataSeries>=levelSignificance[n]) return(paste(rep(userDefSymbol,(length(levelSignificance)-n+1)),sep='',collapse=''))
return('')
}

ListHomogStats<-function(dataYearSeries=NA,refYearSeries=NA,diffFlag=TRUE,homogenization=SNHTabsolute,levelSignificance=c(99,95),criticalValues=climtrends::SNHT.Critical.Values,userDefSymbol='<',nYears=20)
{ # applies an homogenization method and returns a list with the years of a probable shift and other useful results
# equivalent to the "Homogenization Overview" in AnClim
if (any(is.na(dataYearSeries))) stop('<<<dataYearSeries>>> must be a table or dataframe with years and values, no NAs')
if (any(is.na(refYearSeries))) refYearSeries <- NA
if (!is.na(refYearSeries)) if (diffFlag) dataYearSeries[,2] <- dataYearSeries[,2]-refYearSeries[,2] else dataYearSeries[,2] <- dataYearSeries[,2]/refYearSeries[,2]
if (is.function(homogenization)) mystrfun<-as.character(substitute(homogenization)) else stop('<<<homogenization>>> must be a valid function name')
if (mystrfun=='homogenization') {
    nn <- substitute(homogenization)
    i <- 1
    while(TRUE) {
      on <- do.call("substitute", list(as.name(nn), parent.frame(i)))
      if (on == nn) break;
      nn <- on
      i <- i + 1
    }
    mystrfun<-nn
  }
N <- dim(dataYearSeries)[1]
N2 <- dim(dataYearSeries)[2]
Ncrit <- length(levelSignificance)
s1 <- matrix(0, N-1,N2-1)
cat('Series:  Year of change, To value,  Adjust value\n')
for (n in 2:N2) { #n=3
stmp <- GetInfoHomogenization(dataYearSeries=dataYearSeries[,c(1,n)],refYearSeries=NA,nYears=nYears,diffFlag=TRUE,returnData=TRUE,homogenization=homogenization,levelSignificance=levelSignificance,criticalValues=criticalValues)
y2 <- ''
if (!is.null(stmp$endyear)) y2 <- stmp$endyear
cat(n-1,stmp$yearChange, LevelSignifMarks(stmp$ToValue, stmp$criticalFound, userDefSymbol),y2, stmp$ToValue,stmp$AdjustValue,'\n')
# check for breaks at both sides of the "found" break
Ny <- which(dataYearSeries[,1]==stmp$yearChange)
if (Ny>5){
stmpL <- GetInfoHomogenization(dataYearSeries=dataYearSeries[1:(Ny-1),c(1,n)],refYearSeries=NA,nYears=nYears,diffFlag=TRUE,returnData=TRUE,homogenization=homogenization,levelSignificance=levelSignificance,criticalValues=criticalValues)
y2 <- ''
if (!is.null(stmpL$endyear)) y2 <- stmpL$endyear
if (!is.infinite(stmpL$ToValue)) if(stmpL$ToValue> min(stmpL$criticalFound)) cat('\t',n-1,stmpL$yearChange, LevelSignifMarks(stmpL$ToValue, stmpL$criticalFound, userDefSymbol),y2,stmpL$ToValue,stmpL$AdjustValue,
'(',dataYearSeries[Ny,1],'-',dataYearSeries[N,1],',','n =',N-Ny+1,')','\n')
}
if (N-Ny>5){
stmpR <- GetInfoHomogenization(dataYearSeries=dataYearSeries[(Ny):N,c(1,n)],refYearSeries=NA,nYears=nYears,diffFlag=TRUE,returnData=TRUE,homogenization=homogenization,levelSignificance=levelSignificance,criticalValues=criticalValues)
y2 <- ''
if (!is.null(stmpR$endyear)) y2 <- stmpR$endyear
if (!is.infinite(stmpR$ToValue)) if(stmpR$ToValue>min(stmpR$criticalFound)) cat('\t',n-1,stmpR$yearChange, LevelSignifMarks(stmpR$ToValue, stmpR$criticalFound, userDefSymbol),y2,stmpR$ToValue,stmpR$AdjustValue,
'(',dataYearSeries[Ny,1],'-',dataYearSeries[N,1],',','n =',N-Ny+1,')','\n')
}
}
cat('( ',dataYearSeries[1,1],' - ',dataYearSeries[N,1],' ;  n = ',N,' );\t',paste('("',userDefSymbol,'" means that the value exceeds ',min(levelSignificance),'%',ifelse(length(levelSignificance)>1,paste(', "',userDefSymbol,userDefSymbol,'" means ',sort(levelSignificance)[2],'%',sep=''),''),')',sep=''),'(Adjust:from 20 values around the change)\n')
}

ListHomogTZ<-function(dataYearSeries=NA,refYearSeries=NA,diffFlag=TRUE,homogenization=SNHTabsolute, testName='')
{ # List of T statistic and Z statistic for all years
# equivalent to the "Show Details" in AnClim
if (any(is.na(dataYearSeries))) stop('<<<dataYearSeries>>> must be a table or dataframe with years and values, no NAs')
if (any(is.na(refYearSeries))) refYearSeries <- NA
if (!is.na(refYearSeries)) if (diffFlag) dataYearSeries[,2] <- dataYearSeries[,2]-refYearSeries[,2] else dataYearSeries[,2] <- dataYearSeries[,2]/refYearSeries[,2]
N <- dim(dataYearSeries)[1]
Zstat <- scale(dataYearSeries[,2]) # Zi score
s1 <- homogenization(dataYearSeries[,2])
if (is.function(homogenization)) mystrfun<-as.character(substitute(homogenization)) else stop('<<<homogenization>>> must be a valid function name')
if (testName==''){
if (mystrfun %in% c('SNHTabsolute','SNHTrelative')) cat('Single shift of the mean level - Alexandersson Test\nStatistic Ti, Zi and Values\n')
if (mystrfun %in% c('SNHTabsoluteII','SNHTabsoluteII')) cat('Single shift of the mean level II - Alexandersson Test\nStatistic Ti, Zi and Values\n')
if (mystrfun %in% c('SNHTabsolute2SDs','SNHTabsolute2SDs')) cat('Single shift of the mean level and variance - Alexandersson Test\nStatistic Ti, Zi and Values\n')
if (mystrfun %in% c('SNHTabsoluteDoubleShift','SNHTabsoluteDoubleShift')) cat('Double shift of the mean level - Alexandersson Test\nStatistic Ti, Zi and Values\n')
if (mystrfun %in% c('SNHTabsoluteTrend','SNHTabsoluteTrend')) cat('Trend in the mean level - Alexandersson Test\nStatistic Ti, Zi and Values\n')
} else cat(testName,'\n')
if (mystrfun %in% c('SNHTabsolute2SDs','SNHTabsolute2SDs')) return(cbind(year=dataYearSeries[2:(N-2),1],Ti=s1,Zi=Zstat[2:(N-2)],value=dataYearSeries[2:(N-2),2])) else {
if ((mystrfun %in% c('SNHTabsoluteDoubleShift','SNHTabsoluteDoubleShift')) | (mystrfun=='SNHTabsoluteTrend')) return(cbind(year=dataYearSeries[-N,1],Ti=s1$SNHTdata,Zi=Zstat[-N],value=dataYearSeries[-N,2])) else {
if (length(s1)==N) return(cbind(year=dataYearSeries[,1],Ti=s1,Zi=Zstat,value=dataYearSeries[,2])) else return(cbind(year=dataYearSeries[-N,1],Ti=s1,Zi=Zstat[-N],value=dataYearSeries[-N,2]))
}
}
}

"%w/o%" <- function(x,y) x[!x %in% y] #--  x without y, *taken from the R documentation*

FirstMonthOfSeason<-function(intSeason) if (intSeason==1) return (12) else return ((intSeason-1)*3)# returns the first month number for a given season 1 12, 2 3, 36, 4 9
LastMonthOfSeason<-function(intSeason) if (intSeason==1) return (2) else return ((intSeason-1)*3+2)# returns the last month number for a given season 1 2, 2 5, 3 8, 4 11

IsLeapYear<-function(intYear) (intYear %% 400 == 0) | ( (intYear %% 4 == 0) & (intYear %% 100 != 0) )
LastDayOfTheMonth<-function(intMonth,intYear) ifelse(intMonth!=2,c(31,28,31,30,31,30,31,31,30,31,30,31)[intMonth],28+IsLeapYear(intYear))

peaksP<-function(dataSeries,span=3)
{ # autor Brian Ripley
# span has to be odd number
z <- embed(dataSeries, span)
s <- span%/%2
v<- max.col(z) == 1 + s
result <- c(rep(FALSE,s),v)
result <- result[1:(length(result)-s)]
result
}

MonthlyFuncFromDay<-function(yearDF, datecol=1,valcol=2, mfunc=mean)
{
#if (!is.character(yearDF[,datecol])) yearDF[,datecol]<-as.character(yearDF[,datecol])
zu<-unique(format(yearDF[,datecol], "%Y-%m") ) # zu are the unique years
zu.len<-length(zu)
zvec<-vector("numeric",zu.len) # where to store all the years/months
for (n in 1:zu.len) # loop through all years-months
{
m<-1:31
monthdaterange<-as.Date(paste(zu[n],'-',m,sep=''), format="%Y-%m-%d")
zvec[n]<-mfunc(yearDF[which(as.Date(yearDF[,datecol], format="%Y-%m") %in% monthdaterange),valcol], na.rm = T)
}
cbind(year=as.numeric(substr(zu,1,4)),month=as.numeric(substr(zu,6,7)),data=zvec)
}

YearMeanFromDay<-function(yearDF, datecol=1,valcol=2)
{
zu<-unique(format(yearDF[,datecol], "%Y") ) # zu are the unique years
zu.len<-length(zu)
zvec<-vector("numeric",zu.len) # where to store the mean of all years
for (n in 1:zu.len) # loop through all years
# store the mean
zvec[n]<-mean(yearDF[which(yearDF[,datecol] %in% as.Date(paste(zu[n],'-01-01',sep='')):as.Date(paste(zu[n],'-12-31',sep=''))),valcol], na.rm = T)
cbind(as.numeric(zu),zvec)
}

YearFuncFromDay<-function(yearDF, datecol=1,valcol=2, yfunc=mean)
{
#if (!is.character(yearDF[,datecol])) yearDF[,datecol]<-as.character(yearDF[,datecol])
zu<-unique(format(yearDF[,datecol], "%Y") ) # zu are the unique years
zu.len<-length(zu)
zvec<-vector("numeric",zu.len) # where to store the mean of all years
for (n in 1:zu.len) # loop through all years
# store the mean
zvec[n]<-yfunc(yearDF[which(yearDF[,datecol] %in% as.Date(paste(zu[n],'-01-01',sep='')):as.Date(paste(zu[n],'-12-31',sep=''))),valcol], na.rm = T)
cbind(as.numeric(zu),zvec)
}

MonthTrendYearsFuncFromDay<-function(yearDF, datecol=1,valcol=2,mfunc=mean,mmonth)
{
zu<-unique(format(yearDF[,datecol], "%Y") ) # zu are the unique years
zu.len<-length(zu)
zvec<-vector("numeric",zu.len) # where to store all the years/months
for (n in 1:zu.len) # loop through all years-months
{
m<-1:31
monthdaterange<-as.Date(paste(zu[n],'-',mmonth,'-',m,sep=''), format="%Y-%m-%d")
zvec[n]<-mfunc(yearDF[which(as.Date(yearDF[,datecol], format="%Y-%m") %in% monthdaterange),valcol], na.rm = T)
}
cbind(as.Date(paste(zu,'-',mmonth,'-1',sep=''), format="%Y-%m-%d"),zvec)
}

WetDayCount<-function(yearDF, datecol=1,valcol=2,vthreshold)
{
#returns the number of days with precipitation over the threshold
zu<-unique(format(yearDF[,datecol], "%Y") ) # zu are the unique years
zu.len<-length(zu)
zvec<-vector("numeric",zu.len) # where to store the mean of all years
for (n in 1:zu.len) # loop through all years
{
# store the number of days with precipitation over the threshold
tmp<-yearDF[which(yearDF[,datecol] %in% as.Date(paste(zu[n],'-01-01',sep='')):as.Date(paste(zu[n],'-12-31',sep=''))),valcol]
tmp<-tmp[!is.na(tmp)]
zvec[n]<-sum(tmp>=vthreshold)
}
qq<-cbind(as.numeric(zu),zvec)
return(qq)
}

VonNeumannRatio<-function(dataSeries)
{# Von Neumann's Ratio Test for Randomness
n<-length(dataSeries)
m<-mean(dataSeries)
t1<-dataSeries[1:(n-1)]
t2<-dataSeries[2:n]
if (sum((dataSeries-m)^2)==0) return(9999) # to avoid a NaN for a constant series
N<-(sum((t1-t2)^2))/sum((dataSeries-m)^2)
if (is.na(N)) return(9999)
N
}

VonNeumannRatioRank<-function(dataSeries)
{# The Rank Version of Von Neumann's Ratio Test for Randomness
n<-length(dataSeries)
R<-rank(dataSeries)
m<-mean(R)
t1<-R[1:(n-1)]
t2<-R[2:n]
if (sum((dataSeries-m)^2)==0) return(9999) # to avoid a NaN for a division by zero
N<-(sum((t1-t2)^2))/sum((R-m)^2)
if (is.na(N)) return(9999) # to avoid a NaN
N
}

NegCurve<-function(values) 
{#returns -1 for negative shift and 1 for positive shift
if (-max(abs(values), na.rm = T)==min(values, na.rm = T)) 
return(-1) else return(1)
}

GetShiftValue<-function(values)
{#returns the value of the peak (max or min) in a time series
if (NegCurve(values)==1) return(max(values, na.rm = T))
else return(min(values, na.rm = T))
}

ReadYearAnd12MonthsData<-function(filename)
{#reads a monthly values file with the data as y + 12 months in 1 line:
# 1857   -1.7    1.1    4.8   12.3   14.4   16.4   19.8   19.0   14.3   12.9    4.9    1.2
fm1<-read.delim(filename,header=F,sep='',stringsAsFactors=F,colClasses="numeric",na.strings ="-999.99")
z2<-matrix(unlist(fm1[-1]),ncol=13,byrow=F)
z2<-c(apply(z2,1,c))
z3<-rep(fm1[,1],each = 13)
z<-cbind(z3, z2)
z
}

ReadHERTTAdailyCSV<-function(filename)
{#reads a monthly values file from HERTTA daily CSV file
#01.07.1943;0,08;
col.nm<-c('daymonthyear','HERTTAdata','empty') 
z<-read.delim(filename,skip = 4,col.names=col.nm,header =F,stringsAsFactors=F,sep=';',dec = ',')
z<-cbind.data.frame(z,0,0)
z<-z[,-5:-3]
z[,1]<-as.Date(z[,1], "%d.%m.%Y")
z
}

ReadUSHCN<-function(filename)
{#reads a monthly values file from the United States Historical Climatology Network (USHCN)
#01108431895   469    426    571    638    710    789    813    810    786    613  -9999  -9999  -9999 
#STATION ID     1-6   Character
#ELEMENT        7-7   Integer
#YEAR          8-11   Integer
#VALUE1       13-17   Integer
#FLAG1        18-18   Character
#VALUE2       20-24   Integer
#FLAG2        25-25   Character
#  .           .          .
#  .           .          .
#  .           .          .
#VALUE13     97-101   Integer
#FLAG13     102-102   Character
col.widths<-c(6,1,4,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1)
col.names=c('STATION ID','ELEMENT','YEAR','VALUE1','FLAG1','VALUE2','FLAG2','VALUE3','FLAG3','VALUE4','FLAG4','VALUE5','FLAG5','VALUE6','FLAG6','VALUE7','FLAG7','VALUE8','FLAG8','VALUE9
','FLAG9','VALUE10','FLAG10','VALUE11','FLAG11','VALUE12','FLAG12','VALUE13','FLAG13')
zx<-read.fwf(filename,header =F,stringsAsFactors=F,widths=col.widths,na.strings ="-9999",col.names=col.names)
zx}

ReadPSMSLmonthly<-function(filename)
{#reads a monthly values file from the Permanent Service for Mean Sea Level (PSMSL)
#  1943.2084;  6783; 8;000
col.nm<-c('yearmonth','MSL','MissingDays','flag') 
z<-read.delim(filename,col.names=col.nm,header =F,stringsAsFactors=F,sep=';')
z<-cbind.data.frame(z,0,0)
z[,5]<-trunc(z[,1])
z[,6]<-round((z[,1]-z[,5])*12+.5,0)
z[,1]<-as.Date(paste(z[,5],z[,6],1,sep='-'))
z<-z[,c(-5,-6)]
z
}

ReadMeteoSwiss<-function(filename)
{#reads a daily values file from MeteoSwiss
#Year  Month        Temperature      Precipitation  
#1864      1               -5.5               20.6  
col.widths<-c(4,7,19,19)
col.nm<-c('Year','Month','Temperature','Precipitation') 
z<-read.fwf(filename,header =F,stringsAsFactors=F,widths=col.widths,skip = 27,col.names=col.nm)
z[,1]<-gsub(" ", "", paste(z[,1],z[,2],1,sep='-') )
z<-z[,c(-2,-4)]
z[,1]<-as.Date(z[,1],format="%Y-%m-%d")
z[,2]<-as.numeric(z[,2])
return(z) # return the data
}

ReadClimexpKnmiNL12month<-function(filename)
{#reads a monthly values file in climexp.knmi.nl format, a few comments with #, then the data as y + 12 months in 1 line:
# 1857   -1.7    1.1    4.8   12.3   14.4   16.4   19.8   19.0   14.3   12.9    4.9    1.2
# 1829  -11.3  -14.0   -8.9   -1.7    7.5   14.1   17.7   14.5   12.1    4.5   -3.5   -6.6
# 2010  -12.4   -9.1    3.7    4.7 -999.9 -999.9 -999.9 -999.9 -999.9 -999.9 -999.9 -999.9
col.widths<-c(5,7,7,7,7,7,7,7,7,7,7,7,7)
z<-read.fwf(filename,header =F,stringsAsFactors=F,widths=col.widths,comment.char = "#",na.strings ="-999.9")
z2<-matrix(unlist(c(z[,-1])),ncol=12,byrow=F)
z2<-c(apply(z2,1,c))
z3<-rep(z[,1],each = 12)
z<-cbind(z3, z2)
return(z) # return the data
}

ReadClimexpKnmiNL<-function(filename)
{#reads a daily values file in climexp.knmi.nl format, a few comments with #, then the data as y m d:
# 1938  7  1     21.00
col.widths<-c(5,3,3,10)
z<-read.fwf(filename,header =F,stringsAsFactors=F,widths=col.widths,comment.char = "#",na.strings ="-99.99")
#z<-z[!apply(is.na(z), 1, any),]# remove all rows with at least 1 NA on a column
z[,1]<-as.Date(paste(z[,1],z[,2],z[,3],sep='-'),format="%Y-%m-%d")
z<-z[,c(-2,-3)]
return(z) # return the data
}

ReadCRDEC<-function(filename, FlagAnnualMean=F)
{
#reads a file from GHCN-D, with date+Q into a dataframe
#61333601895 40-99999M-99999M-99999M-99999M-99999M-99999M-99999M-99999M-99999M-99999M-99999M-99999M
col.widths<-c(7,4,3,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1,6,1)
z<-read.fwf(filename,header =F,stringsAsFactors=F,na.strings ="-99999",widths=col.widths)
if (!FlagAnnualMean) return (z[-c(1,5,7,9,11,13,15,17,19,21,23,25,27)])
else
{
z<-z[-c(1,5,7,9,11,13,15,17,19,21,23,25,27)]
z[,2]<-apply(z[,3:14],1,mean)
return(z[-(3:14)])
}
}

ReadGHCNymd<-function(filename)
{#Global Historical Climatology Network-Daily (GHCN-Daily)
#reads a file from GHCN-D, with date+Q into a dataframe
# 1917  4  6      5.00
col.widths<-c(5,3,3,10)
z<-read.fwf(filename,header =F,stringsAsFactors=F,skip = 4,widths=col.widths)
z<-z[!apply(is.na(z), 1, any),]# remove all rows with at least 1 NA on a column
z[,1]<-as.Date(paste(z[,1],z[,2],z[,3],sep='-'),format="%Y-%m-%d")
z<-z[,c(-2,-3)]
return(z) # return the data
}

ReadHERTTAdmyQdaily<-function(filename)
{#reads a file from HERTAA *Finland*, with date+Q into a dataframe
z<-read.delim(filename,header =F,stringsAsFactors=F)
z[,1]<-as.Date(as.character(z[,1]), format="%d.%m.%Y")#turn it into a date datatype
return(z) # return the data
}

ReadEtmgegFile<-function(filename,optCols=NULL)
{#reads a Royal Netherlands Meteorological Institute (KNMI) file into a dataframe
col.nm<-c('STN','YYYYMMDD','DDVEC','FHVEC','FG','FHX','FHXH','FHN','FHNH','FXX','FXXH','TG','TN','TNH','TX','TXH','T10N','T10NH','SQ','SP','Q','DR','RH','RHX','RHXH','PG','PX','PXH','PN','PNH','VVN','VVNH','VVX','VVXH','NG','UG','UX','UXH','UN','UNH','EV24','junk')
z<-read.csv(filename,header =F,skip = 49,col.names=col.nm)
if (length(optCols)>0) z<-z[sort(optCols)] #user select colums
z<-z[!apply(is.na(z), 1, any),]# remove all rows with at least 1 NA on a column
#if the user selected the date column (2), turn it into a date datatype
if (length(optCols)>0) if (2 %in% optCols) 
if (1 %in% optCols) z[,2]<-as.Date(as.character(z[,2]), format="%Y%m%d")
else z[,1]<-as.Date(as.character(z[,1]), format="%Y%m%d")
return(z) # return the data
}

ReadUSGS2timeseries <-function(filename)
{# reads a .part file which is a USGS discharge file downloaded from 
#example: t<-ReadUSGS2timeseries("03275000.part")
z<-read.table(filename,header = T,comment.char="#",stringsAsFactors=F)
z<-z[-1,] # remove 1st line, it has no data
z<-z[,-1:-4] # keep only the last 2 columns, year and discharge - year_nu mean_va
z<-as.matrix(cbind(as.numeric(z[,1]),as.numeric(z[,2]))) # convert to a numeric matrix
ts.y<-ts(z[,2],start = min(z[,1]), frequency = 1) # convert to a time series
}

AllReadUSGS2timeseries <-function()
{# reads multiple USGS discharge files
lf<-list.files(path = ".", pattern ="\\.part")
mylist<-vector("list", length(lf))
names(mylist)=lf
for (lfile in lf)
mylist[lfile ]<-list(ReadUSGS2timeseries(lfile))
return(mylist)
}

ReadECAdata<-function(filename)
{#European Climate Assessment & Dataset project
#FILE FORMAT (MISSING VALUE CODE IS -9999):
#01-06 SOUID: Source identifier
#08-15 DATE : Date YYYYMMDD
#17-21 RR   : daily precipitation amount in 0.1 mm
#23-27 Q_RR : quality code for RR (0='valid'; 1='suspect'; 9='missing')
z<-read.csv(filename,header =F,skip = 21,col.names=c('SOUID','DATE','RR','Q_RR'))
z<-z[,-1]# remove the SOUID
z<-z[which(z$Q_RR==0),1:2]#keep valid data only, 0='valid'
z[,1]<-as.Date(paste(substr(z[,1],1,4),substr(z[,1],5,6),substr(z[,1],7,8),sep='/'), format="%Y/%m/%d")
z
}

ReadECAdataIndexDTR<-function(filename)
{#European Climate Assessment & Dataset project
# DATA IS FOR INDEX Mean of diurnal temperature range (DTR) with unit 0.01 Temperature deg. C
# FILE FORMAT (MISSING VALUE CODE = -999999):
# 01-06 SOUID: Source identifier
# 08-11 YEAR : YYYY
# 13-20 ANNUAL DATA VALUES
# 21-28 WINTER HALF YEAR DATA VALUES
# 29-36 SUMMER HALF YEAR DATA VALUES
# 37-44 WINTER (DJF) DATA VALUES
# 45-52 SPRING (MAM) DATA VALUES
# 53-60 SUMMER (JJA) DATA VALUES
# 61-68 AUTUMN (SON) DATA VALUES
# 69-76 JANUARY DATA VALUES
# etc. 
# 157-164 DECEMBER DATA VALUES
z<-read.delim(filename,sep='',header =F,skip = 30,col.names=c('SOUID','YEAR','ANNUAL DATA','WINTER HALF YEAR DATA','SUMMER HALF YEAR DATA',
'WINTER (DJF) DATA','SPRING (MAM) DATA','SUMMER (JJA) DATA','AUTUMN (SON) DATA','JANUARY','FEBRUARY','MARCH','APRIL','MAY','JUNE',
'JULY','AUGUST','SEPTEMBER','OCTOBER','NOVEMBER','DECEMBER'))
z<-z[,-1]# remove the SOUID
z
}

ReadNdp040stationInventory<-function(filename)
{
#reads the contents of station_inventory.txt into a table
#example: ndpStationInventory<-ReadNdp040stationInventory('../Russian stations/ndp040/station_inventory.txt')
col.widths<-c(6,25,7,7,8,5,5,5,5,5,5,5,5,4)
col.nm<-c('WMO', 'NAME', 'LAT', 'LON', 'ELEV', 'MINTFYR', 'MINTMISS', 'MIDTFYR', 'MIDTMISS', 'MAXTFYR', 'MAXTMISS', 'PRCPFYR', 'PRCPMISS', 'LYR')
z<-read.fwf(filename,header =F,stringsAsFactors=F,col.names=col.nm,widths=col.widths)
return(z)
}

ReadNdp040stationHistory<-function(filename)
{
#reads the contents of station_history.txt into a table
#example: ndpStationHistory<-ReadNdp040stationHistory('../Russian stations/ndp040/station_history.txt')
col.widths<-c(5,5,5,3,3,3,3)
col.nm<-c('WMO', 'TYPE', 'YEAR', 'MONTH', 'DAY', 'DIST', 'DIRECT')
z<-read.fwf(filename,header =F,stringsAsFactors=F,col.names=col.nm,widths=col.widths)
z$YEAR[which(z$YEAR==-999)]<-NA
z$MONTH[which(z$MONTH==-9)]<-NA
z$DAY[which(z$DAY==-9)]<-NA
return(z)
}

ReadNdp040stationDat<-function(filename)
{
#reads the contents of a station data file (.dat) into a table
#example: ndpStationDat<-ReadNdp040stationDat('../Russian stations/ndp040/f20674.dat')
col.widths<-c(6,6,4,4,2,6,2,6,2,6,2,6,2,3,4)
col.nm<-c('WMO', 'YEAR', 'MONTH', 'DAY','TFLAG', 'TMIN', 'QTMIN', 'TMID', 'QTMID', 'TMAX', 'QTMAX', 'R', 'CR', 'QR', 'DATA_FLAGS')
z<-read.fwf(filename,header =F,stringsAsFactors=F,col.names=col.nm,widths=col.widths)
return(z)
}

ReadNdp040stationData<-function(filename)
{
#reads the contents of a station data file (.data) into a table
#example: ndpStationData<-ReadNdp040stationData('../Russian stations/ndp040/ussr1.data')
col.widths<-c(5,4,4,2,2,2,4,1,1,2,4,1,1,2,4,1,1)
col.nm<-c('WMO', 'TYPE', 'YEAR', 'MONTH', 'NOBS', 'DAY4', 'DATA4', 'FLAGA4', 'FLAGB4', 'DAY6', 'DATA6', 'FLAGA6', 'FLAGB6', 'DAY15', 'DATA15', 'FLAGA15', 'FLAGB15'
)
z<-read.fwf(filename,header =F,stringsAsFactors=F,col.names=col.nm,widths=col.widths)
return(z)
}

GetNDP<-function(ndpDat, climVar)
{
# returns one climate variable from one NDP file read into a data.frame
# climVar = TMIN, TMID, TMAX, R, DTR R is daily precipitation total
if (climVar %in% c('TMIN', 'TMID', 'TMAX', 'R')) z<-ndpDat[,c('YEAR', 'MONTH', 'DAY', climVar, paste('Q',climVar,sep=''))]
if (climVar =='DTR') z<-ndpDat[,c('YEAR', 'MONTH', 'DAY', 'TMIN', 'QTMIN', 'TMAX', 'QTMAX')]
if (climVar=='TMID') z<-z[which(z$QTMID==0),] # keep only valid values 
if (climVar=='TMIN') z<-z[which(z$QTMIN==0),] # keep only valid values 
if (climVar=='TMAX') z<-z[which(z$QTMAX==0),] # keep only valid values 
if (climVar=='R') z<-z[which(z$QR==0),] # keep only valid values 
if (climVar=='DTR') z<-z[which(z$QTMAX==0 & z$QTMIN==0),] # keep only valid values 
z<-z[-5] # remove Q flag
if (climVar=='DTR')
{
z<-z[-6] # remove Q flag
z<-cbind(z,abs(z[5]-z[4])) # DTR=MAX-MIN
z<-z[-4];z<-z[-4] # remove MAX MIN
}
z[,1]<-as.Date(as.character(paste(z[,3],z[,2],z[,1],sep='-')), format="%d-%m-%Y")#turn it into a date datatype
z<-z[c(-2,-3)] # remove month and day data
z
}

AllfilesIEObyPathPattern <-function(mypath, mypattern)
{#reads a daily values file from IEO (Instituto Español de Oceanografía, Spain)
#  2Santande   LAT=43 28  N LONG=003 48  W   REF: TGZ
#  year   month   day    hour    heigh(m)   station
#  1944      1      1      0       1.080       2
flg<-T
lf<-list.files(path = mypath, pattern = mypattern)
lf
z2<-0
for (lfile in lf){
#print(lfile)
col.widths<-c(6,7,7,7,12,8)
col.nm<-c('year','month','day','hour','heigh','station') 
z<-read.fwf(paste(mypath,lfile ,sep='/'),widths=col.widths,col.names=col.nm,stringsAsFactors=F, skip = 2)
if (flg) z2<-z
else z2<-rbind.data.frame(z2,z)
flg<-F
}
z2[,1]<-gsub(" ", "", paste(z2[,1],z2[,2],z2[,3],sep='-') )
z2[,1]<-as.Date(z2[,1],format="%Y-%m-%d")
z2<-z2[,c(-2,-3,-4,-6)]
return(z2)
}

SiegelTukey<-function(x,y,id.col=FALSE,adjust.median=FALSE,rnd=-1,alternative="two.sided",mu=0
,paired=FALSE,exact=FALSE,correct=TRUE,conf.int=FALSE,conf.level=0.95,showresult=TRUE, returnresult=TRUE){
# 
  if(id.col==FALSE){
    data=data.frame(c(x,y),rep(c(1,2),c(length(x),length(y))))
    } else {
    data=data.frame(x,y)
    }
  names(data)=c("x","y")
  data=data[order(data$x),]
  if(rnd>-1){data$x=round(data$x,rnd)}
  if(adjust.median==T){
 data$x[data$y==1]=data$x[data$y==1]-(median(data$x[data$y==1])-median(data$x[data$y==2]))/2
 data$x[data$y==2]=data$x[data$y==2]-(median(data$x[data$y==2])-median(data$x[data$y==1]))/2
  }
  md1 <- median(data$x[data$y==1])
  if (showresult) cat("Median of group 1 = ",md1,"\n")
  md2 <- median(data$x[data$y==2])
  if (showresult) cat("Median of group 2 = ",md2,"\n","\n")
  tmd <- wilcox.test(data$x[data$y==1],data$x[data$y==y])
  if (showresult) cat("Test of median differences\n",unlist(wilcox.test(data$x[data$y==1],data$x[data$y==y])),"\n")
  a=rep(seq(ceiling(length(data$x)/4)),each=2)
  b=rep(c(0,1),ceiling(length(data$x)/4))
  rk.up=c(1,(a*4+b))[1:ceiling(length(data$x)/2)]
  rk.down=rev(c(a*4+b-2)[1:floor(length(data$x)/2)])
  if (showresult) cat("Performing Siegel-Tukey rank transformation...","\n","\n")
  rks=c(rk.up,rk.down)
  unqs=unique(sort(data$x))
  corr.rks=tapply(rks,data$x,mean)
  cbind(unqs,corr.rks)
  rks.data=data.frame(unqs,corr.rks)
  names(rks.data)=c("unique values of x","tie-adjusted Siegel-Tukey rank")
  print(rks.data,row.names=F)
  names(rks.data)=c("unqs","corr.rks")
  data=merge(data,rks.data,by.x="x",by.y="unqs")
  rk1=data$corr.rks[data$y==1]
  rk2=data$corr.rks[data$y==2]
  if (showresult) cat("\n","Tie-adjusted Siegel-Tukey ranks of group 1","\n")
  group1=data.frame(data$x[data$y==1],rk1)
  names(group1)=c("x","rank")
  if (showresult) print(group1,row.names=F)
  if (showresult) cat("\n","Tie-adjusted Siegel-Tukey ranks of group 2","\n")
  group2=data.frame(data$x[data$y==2],rk2)
  names(group2)=c("x","rank")
  if (showresult) print(group2,row.names=F)
  if (showresult) cat("\n")
  if (showresult) cat("Siegel-Tukey test","\n")
  if (showresult) cat("Siegel-Tukey rank transformation performed.","Tie adjusted ranks computed.","\n")
  if(adjust.median==T)cat("Medians adjusted to equality.","\n") else
  if (showresult) cat("Medians not adjusted.","\n")
  if (showresult) cat("Rank sum of group 1 =", sum(rk1),"    Rank sum of group 2=",sum(rk2),"\n")
  w <- wilcox.test(rk1,rk2,alternative=alternative,mu=mu,paired=paired,exact=exact,correct=correct,conf.int=conf.int,conf.level=conf.level)
 if (showresult) print(wilcox.test(rk1,rk2,alternative=alternative,mu=mu,paired=paired,exact=exact,correct=correct,conf.int=conf.int,conf.level=conf.level))
if (returnresult) return(list(groupmedian1=md1,groupmedian2=md2,testmediandifferences=tmd,wilcox.test.p.value=tmd$p.value,
wilcox.test.null.value=tmd$null.value,unique.x.tieadjusted.rank=rks.data,tieadjusted.ranks.group1=group1,tieadjusted.ranks.group2=group2,wilcoxon.W=w[["statistic"]],wilcoxon.p.value=w$p.value))
 }

VDTR<-function(yearDF,datecol=1,valcol=2)
{ # annual mean of the absolute day-to-day differences of the diurnal temperature range (VDTR) 
zu<-unique(format(yearDF[,datecol], "%Y") ) # zu are the unique years
zu.len<-length(zu)
zvec<-vector("numeric",zu.len) # where to store the mean of all years
for (n in 1:zu.len) # loop through all years
# store the mean of the absolute day-to-day differences of the DTR
{
b<-0
ytmp<-yearDF[which(yearDF[,datecol] %in% as.Date(paste(zu[n],'-01-01',sep='')):as.Date(paste(zu[n],'-12-31',sep=''))),valcol]
for (k in 2:length(ytmp))
b<-b+abs(ytmp[k]-ytmp[k-1])
b<-b/(length(ytmp)-1)
zvec[n]<-b
}
zvec
}

SNHTabsolute<-function(dataSeries)
{#calculates the SNHT absolute, assumes no NAs
tsn<-(dataSeries-mean(dataSeries))/sd(dataSeries)
n<-length(dataSeries)
zvec<-vector("numeric",n-1) # to store the result
for (v in 1:(n-1)) zvec[v]<-v*(mean(tsn[1:v]))^2+(n-v)*(mean(tsn[(v+1):n]))^2
zvec
}

SNHTabsoluteII<-function(dataSeries)
{#calculates the SNHT absolute, standard deviation different from 1,
# assumes no NAs, (A4) and (A5) from HOMOGENIZATION OF SWEDISH TEMPERATURE DATA. PART I: HOMOGENEITY TEST FOR LINEAR TRENDS
tsn<-(dataSeries-mean(dataSeries))/sd(dataSeries)
n<-length(dataSeries)
zvec<-vector("numeric",n-1) # to store the result
for (v in 1:(n-1)) zvec[v]<- -2*n*log( sqrt(( n-1-( v*(mean(tsn[1:v]))^2+(n-v)*(mean(tsn[(v+1):n]))^2 ) )/n) ) -1
zvec
}

SNHTabsolute2SDs<-function(dataSeries)
{#calculates the SNHT absolute, two different standard deviations, assumes no NAs
# HOMOGENIZATION OF SWEDISH TEMPERATURE DATA. PART I: HOMOGENEITY TEST FOR LINEAR TRENDS
tsn<-(dataSeries-mean(dataSeries))/sd(dataSeries)
#tsn<-dataSeries
n<-length(dataSeries)
zvec<-vector("numeric",n-3) # to store the result
for (v in 2:(n-2)) zvec[v-1]<- -2*v*  log(sqrt(( sum(tsn[1:v]^2) - sum(tsn[1:v])^2/v)/v)) - 2*(n-v)*  log(sqrt(( sum(tsn[(v+1):n]^2) - sum(tsn[(v+1):n])^2/(n-v))/(n-v)))-1
zvec
}

SNHTabsoluteDoubleShift<-function(dataSeries)
{# Double shift of the mean level
tsn<-(dataSeries-mean(dataSeries))/sd(dataSeries)
n<-length(dataSeries)
MatS <- matrix(-Inf,(n-1),n) # to store the result
for (a in 1:(n-1)) for (b in (a+1):n) {
MatS[a,b] <- a*(mean(tsn[1:a]))^2 + (b-a)*(mean(tsn[(a+1):b]))^2 + (n-b)*(mean(tsn[( ifelse(b==n,0,b) +1):n]))^2
}
SNHTret <- apply(MatS,1,max)
list(SNHTdata=SNHTret, PosSecondShift=apply(MatS,1,function(x) which(x == max(x))))
}

SNHTabsoluteTrend<-function(dataSeries)
{#calculates the SNHT absolute with trend in the mean level, assumes no NAs
tsn<-(dataSeries-mean(dataSeries))/sd(dataSeries)
n<-length(dataSeries)
MatS <- matrix(-Inf,(n-1),n) # to store the result
for (a in 1:(n-1)) for (b in (a+1):n) {
Z1m <- mean(tsn[1:a])
Z2m <- ifelse(b==n,0,mean(tsn[(b+1):n]))
SA=sum(((a+1):b-a)^2/(b-a)^2)
SB=sum((b-(a+1):b)^2/(b-a)^2)
SZA=sum(tsn[(a+1):b]*((a+1):b-a)/(b-a))
SZB=sum(tsn[(a+1):b]*(b-(a+1):b)/(b-a))
SAB=sum((b-(a+1):b)*((a+1):b-a)/(b-a)^2)
SK=-SAB / (SA+n-b)
SL=((n-b)*Z2m+SZA) / (SA+n-b)
u1=(a*Z1m+SZB-SL*SAB) / (a+SB+SK*SAB)
u2=u1*SK+SL
MatS[a,b] <- -a*u1^2 + 2*a*u1*Z1m - u1^2*SB - u2^2*SA + 2*u1*SZB + 2*u2*SZA - 2*u1*u2*SAB - (n-b)*u2^2 + 2*(n-b)*u2*Z2m
}
SNHTret <- apply(MatS,1,max)
list(SNHTdata=SNHTret, PosStartTrend=apply(MatS,1,function(x) which(x == max(x))))
}

BuishandRangeTest<-function(dataSeries)
{
n<-length(dataSeries)
s2<-1:n
s<-sapply(s2,function(k) {
m<-mean(dataSeries)
b<-sum(dataSeries[1:k]-m)
return(b)
})
s/sd(dataSeries)/sqrt(n)
}

PettittTest<-function(dataSeries)
{
n<-length(dataSeries)
s2<-1:n
s<-sapply(s2,function(k) {
n<-length(dataSeries)
r<-rank(dataSeries)
b<-0
for (kk in 1:k)
b<-b+r[kk]
return(2*b-kk*(n+1))
})
s
}

CraddockTest <- function(yearDF, valcol1,valcol2)
{
meanC1 <-mean(yearDF[,valcol1])
meanC2 <-mean(yearDF[,valcol2])
c2C <- yearDF[,valcol2]+meanC1-meanC2
cDiff <- yearDF[,valcol1] - c2C
cumsum(cDiff)
}

WorsleyLikelihoodRatio <- function(yearDF, datecol=1,valcol=2, returnZk=FALSE)
{# Worsley Likelihood Ratio (parametric test for step jump in mean)
# Trend 1.0.2 User Guide, chapter 4.2.6 Worsley Likelihood Ratio Test, pp. 19
Sk <- cumsum(yearDF[,valcol]-mean(yearDF[,valcol]))
x.m2 <- (yearDF[,valcol]-mean(yearDF[,valcol]))^2
n <- dim(yearDF)[datecol]
N <- 1:n
Zk <- Sk/sqrt(N * (n-N))
x.m2.C <- sqrt(sum(x.m2)/n)
Zk <- Zk/x.m2.C
Zk2 <- abs(Zk)
Zk2 <- Zk2[is.finite(Zk2)]
V <- max(Zk2)
YearChange<-(yearDF[which(Zk2==V),datecol])
W <- sqrt(n-2)*V/sqrt(1-V^2)
if (!returnZk) return(list(YearChange=YearChange,V=V,W=W)) else return(list(YearChange=YearChange,V=V,W=W, Zk=Zk))
}

