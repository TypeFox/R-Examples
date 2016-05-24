Gaptest<-function(DesY) {
# function to compute gap statistic
ncheck<-dim(DesY)
ncheck<-ncheck[1]
 tcnd=TRUE
 if (ncheck==8) {tcnd=FALSE}
 if (ncheck==16) {tcnd=FALSE}
 if (ncheck==32) {tcnd=FALSE}
 if (tcnd) {stop("This function only works for 8, 16, or 32 run designs","\n")
    } else {
if (ncheck==8) ncheck=16
#####################################
# 50th and 99th percentiles of the gap statistic ##
critg16<-c(1.7884,5.1009)
critg32<-c(1.7297,5.8758)

### First Pass through the data ####

###### Step 1  #######
#fit model to saturated design
modf<-lm(y~(.)^4,x=TRUE,data=DesY)


#extract the regression coefficients
nbeta<-dim(DesY)
nbeta<-nbeta[1]
he<-modf$coef
# This extracts the coefficients that are not NA
selcol<-which(!is.na(he))
he<-he[selcol]
he<-he[-1]
#number of coefficients
p<-length(he)
#number of runs
n<-p+1
# This trims unnecessary characters from coefficient names
cn1<-names(he)
ccn1<-gsub("[^A-Z]","",cn1)
names(he)<-ccn1

##### End of Step 1 ########

###### Steps 2 and 3  #######
#calculate the pse statistic
ahe<-abs(he)
s0<-1.5*median(ahe)
selhe<-ahe<(2.5*s0)
pse=1.5*median(ahe[selhe])
#library(BsMD)
#pse<-LenthPlot(modf,plt=FALSE)
#pse<-pse[2]
#calculate the gap statistic
gap<-gapstat(he,pse)
# checks to see if gap statistic exceeds 50th percentile
if (ncheck==16) {test=(gap>critg16[1])
 } else  {test=(gap>critg32[1])}
##### End Step 2 and 3  #####

if (test) {
##### Step 4  #####
#extract the model X matrix
X<-modf$x
# This selects columns of the X matrix that correspond to non-missing
# coefficients
X<-X[,selcol]
X<-X[,-1]
#gets signs of regression coefficients
se<-as.matrix(sign(he),nrow=1)
# find signigicant effects using LGB
###### check he
#cat("heP23 = ",he,"\n")
###########
LGB(he)
sigef<-LGBc(he,rpt=FALSE,plt=FALSE)
##########
#cat("sigefP23=",sigef,"\n")
##############
# make signs of significant effects zero
 for (i in 1:length(he)) {
     if (sigef[i]=="yes")   {se[i]=0 }
                         }
#gets sum of products of signed effects and rows of X matrix
sp<-X%*%se

#finds index of largest sum of products as index of potential outlier
asp<-abs(sp)
oo<-max.col(t(asp))

### End Step 4 ####

###### Step 5  ######

# calculates the bias
 # first get absolute regression coefficients
ae<-abs(he)
 # next sort absolute effects
sae<-sort(ae)
 #get the number of effects in smallest half
nsmall<-round(length(he)/2)
 # sum the smallest half absolute effects to get bias
bias<-2*sum(sae[1:nsmall])

##### Step 6 ######
 # gets corrected response vector
y<-DesY$y
ycorr<-DesY$y
ycorr[oo]<-ycorr[oo]+(-1*sign(sp[oo]))*bias
 # makes vector of indicators for outlier
detect<-c(rep("no",n))
detect[oo]<-"yes"
cat("Initial Outlier Report","\n")
cat("Standardized-Gap = ",gap, "Significant at 50th percentile","\n")
### End of first pass throught the data #######

### Second Pass throught the data ###########
### Step 1 ####
# augment DesY with corrected data
DesYc<-cbind(DesY[,1:(dim(DesY)[2]-1)],ycorr)
# fit saturated model to corrected data
modf<-lm(ycorr~(.)^4,x=TRUE,data=DesYc)

#extract the regression coefficients
che<-modf$coef
# This extracts the coefficients that are not NA
che<-che[!is.na(che)]
che<-che[-1]
#number of coefficients
p<-length(che)
#number of runs
n<-p+1
# This trims unnecessary characters from coefficient names
cn<-names(che)
ccn<-gsub("[^A-Z]","",cn)
names(che)<-ccn
### End of Step 1 ####

###### Steps 2 and 3  #######
#calculate the pse statistic
ache<-abs(che)
s0<-1.5*median(ache)
selche<-ache<(2.5*s0)
psec=1.5*median(ache[selche])

#psec<-LenthPlot(modf,plt=FALSE)
#psec<-psec[2]
#calculate the gap statistic
gap<-gapstat(he,psec)
# checks to see if gap statistic exceeds 99th percentile
if (ncheck==16) test2=(gap>critg16[2]) else  test2=(gap>critg32[2])
##### End Step 2 and 3  #####

if (test2) {
cat("Final Outlier Report","\n")
cat("Standardized-Gap = ",gap, "Significant at 99th percentile","\n")
cat("   ","\n")
cat("    Corrrected Data Report  ","\n")
cat("Response  Corrected Response   Detect Outlier","\n")
cat(paste(format(DesY$y, width=8), format(DesYc$ycorr, width=13),"           ", format(detect, width=10),"\n"),sep="")

# use LGBc to test significance of effects calculated from corrected data

tce<-LGBc(che)
  } else {
cat("Final Outlier Report","\n")
cat("No significant outlier detected in second pass","\n" )
# use LGB to test significance of effects calculated from corrected data
LGB(he)
cat("    ","\n")
         }

### End of second pass through the data #####
  }
}
# end of function Gaptest
 }

