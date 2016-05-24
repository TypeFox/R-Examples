plend <-
function(dd){
# DESCRIPTION
#
# Compute Product Limit Estimate(PLE) of F(x) for positive data with 
#  non-detects ( left censored data) 
#
# USAGE
# plend(dd)
#
# ARGUMENTS
# dd is an n by 2 matrix or data frame
# dd[,1]=  exposure variable in column 1
# dd[,2] = ( 0 for non-detect 1 for detect)
#
# DETAILS  
# Direct calculation of PLE for left censored data as defined by
# Schmoyer et al (see ORNL/TM-2005/52 Section 3.5 for details)
#
# VALUE
# a(j)    is the  value of jth detect (ordered)
# ple(j)  is PLE of F(x) at a(j) 
# n(j)    number of detects or non-detects <= a(j)
# r(j)    number of detects equal to a(j)
# surv(j) 1 - ple(j) is PLE of S(x) 
# 
# REFERENCES
# Schmoyer,Beachaump,Brand, and Hoffman et al (1996)"
# 
#              Environmental and Ecological  Statistics,3 81-79
#
# Frome, E. L. and Wambach, P.F. "
#  ORNL/TM-2005/52 Section 3.5
#
# SEE ALSO:
#  pleicf, plekm, bootple 
#
# COMMENTS
#
# In survival analysis S(x)= 1 - F(x) is the survival function
# i.e S(x) = P[X> x]. In  environmental and occupational situations
# 1 - F(x) is the "exceedance" function, i.e C(X) = 1-F(x) = P[X>=x]
#  
# EXAMPLE
# data(aihand) ## use AIHA data set with nondetects
# ple<- plend(aihand)
# qqpln(ple) #  lognormal q-q plot based on PLE for censored data
# 
nn<- length(dd[,1]); n<-rev( 1: nn )
# If xmnd=TRUE and MAXIMUM X IS A NON-DETECT CHANGE IT TO A DETECT
dd<- dd[ order(dd[,1],dd[,2]),]

if( dd[nn,2]==0 ){ print("maximium x is non-detect") }
#         if(xmnd) {dd[nn,2]<-1 ;print("CHANGE IT TO A DETECT")}}

#  dd[,2] in order so that detects preceed non-detects for each a(j)
t1<- as.matrix( cbind( dd[ rev(order(dd[,1],dd[,2])) ,1:2],n ))
#       select detects and calculate n(j) and r(j) for each a(j)
t2<- t1[ t1[,2]==1,]
ung<- c(-1, diff(t2[,1] )) 
d<- rev( table( t2[,1] ) ) 
#  columns of t2 are x[j] n[j] r[j] for det=1 IN REVERSE order
if(length(d)==1) t2<-matrix( c(t2[ ung< 0,] ,d),1 )[c(1,3:4)] 
else t2<- as.matrix( cbind( t2[ ung< 0,],d )[,c(1,3:4)] )
#  x contains detects and nx is number of unique detects
x<- dd[ dd[,2]==1,1] ; nx<-length(x) ; ndt<- max(x)
# add a row for F(ao) if min(x) is a non-detect
if( nx < dim(dd)[1] ){
   ndt<- dd[ dd[,2]==0, 1] #  nondetects
   if( min(ndt) < min(x) ){
   nl<- sum( ifelse(dd[,1] <= min(ndt),1,0) )
   frow<- c(min(ndt),nl,0)
   t2<- rbind(t2,frow)
   }
}
if( dim(t2)[1] < 2){print(t2); stop("At least 2 detects required for PLE")}
#   calculate the ple = Prod( [n[j] - d[j] ]/n[j] ) 
n1<- dim(t2)[1]
ple<- cumprod( (t2[,2]-t2[,3])/t2[,2] )
ple<- c(1,ple[1:n1-1])  # PLE of F(x)
suv<- 1 - ple           # PLE of S(x) (Kaplan-Meier)
t2<- cbind(t2,ple,suv)
t2<- t2[order(t2[,1]),c(1,4,2,3,5)]
dimnames(t2)<-list( NULL,c("a","ple","n", "r","surv") )
t2<-data.frame(t2)
t2
}

