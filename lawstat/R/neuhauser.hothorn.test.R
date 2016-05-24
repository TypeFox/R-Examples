`neuhauser.hothorn.test` <-
function(y, group, location = c("median", "mean", "trim.mean"), tail = c("right","left","both"), 
trim.alpha = 0.25, bootstrap = FALSE, num.bootstrap = 1000, 
correction.method = c("none","correction.factor","zero.removal","zero.correction"))
{

### stop the code if the length of y does not match the length of group ###

 if (length(y)!=length(group))
 {
  stop("the length of the data (y) does not match the length of the group")
 }

### install the mvtnorm package ###

#library(mvtnorm)

### assign location, tail, and a correction method ###

 location <- match.arg(location)
 tail <- match.arg(tail)
 correction.method <- match.arg(correction.method)
 DNAME = deparse(substitute(y))
 y<-y[!is.na(y)]
 group<-group[!is.na(y)]  

### stop the code if the location "trim.mean" is selected and trim.alpha is too large ###

 if ((location == "trim.mean") & (trim.alpha == 1)) 
 {
  stop("trim.alpha value of 0 to 0.5 should be provided for the trim.mean location")
 }
 
### sort the order just in case the input is not sorted by group ###

 reorder <- order(group)
 group <- group[reorder]
 y<-y[reorder]
 gr<-group
 group<-as.factor(group)
 original.group<-group

### define the measure of central tendency (mean, median, trimmed mean) ###

 if (location == "mean") 
 {
  means <- tapply(y, group, mean)
  METHOD="double contrast test based on the absolute deviations from the mean"
 }

 else if (location == "median") 
 {
  means <- tapply(y, group, median)
  METHOD="double contrast test based on the absolute deviations from the median"
 }

 else 
 {
  location = "trim.mean"
  trimmed.mean <- function(y) mean(y, trim = trim.alpha)
  means <- tapply(y, group, trimmed.mean)
  METHOD="double contrast test based on the absolute deviations from the trimmed mean"
 }

### calculate the sample size of each group and absolute deviation from center ###   
 
 n<-tapply(y, group, length)
 original.n<-n
 ngroup<-n[group]
 resp.mean<-abs(y - means[group])

### assign no correction technique if the central tendency is median, and ###
### any technique other than "correction.factor" is chosen                ###

 if(location != "median" && correction.method != "correction.factor")
 { 
  METHOD<-paste(METHOD,"(",correction.method,"not applied because the location is not set to median",")")
  correction.method = "none"
 }

### multiply the correction factor to each observation if "correction.factor" is chosen ###

 if(correction.method == "correction.factor") 
 {
  METHOD<-paste(METHOD,"with correction factor applied")
  correction<-1/sqrt(1-1/ngroup)
  resp.mean<-resp.mean*correction
  r.means <- tapply(resp.mean,group,mean)
 }

### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
### or "zero.correction" (Noguchi and Gel, 2009).                            ###

 else if(correction.method == "zero.correction"||correction.method == "zero.removal")
 {
  if(correction.method == "zero.removal")
  {
   METHOD<-paste(METHOD,"with Hines-Hines structural zero removal method")
  }
  if(correction.method == "zero.correction")
  {
   METHOD<-paste(METHOD,"with modified structural zero removal method and correction factor")
  }

### set up variables for calculating the deviation from center ###

  mean.diff <- y-means[group]
  k<-length(n)
  temp<-double()
  endpos<-double()
  startpos<-double()

### calculate the absolute deviation from mean and remove zeros ###

  for(i in 1:k)
  {
   group.size<-n[i]
   j<-i-1

### calculate the starting and ending index of each group ###

   if(i==1)start<-1

   else start<-sum(n[1:j])+1

   startpos<-c(startpos,start)
   end<-sum(n[1:i])
   endpos<-c(endpos,end)

### extract the deviation from center for the ith group ###

   sub.mean.diff<-mean.diff[start:end]
   sub.mean.diff<-sub.mean.diff[order(sub.mean.diff)]

### remove structural zero for the odd-sized group ###

   if(group.size%%2==1)
   {
    mid<-(group.size+1)/2
    temp2<-sub.mean.diff[-mid]
   }

### remove structural zero for the even-sized group ###

   if(group.size%%2==0)
   {
    mid<-group.size/2
    mid2<-mid+1

### set up the denominator value for the transformation ###

### set 1 for the "zero.correction" option ###

    denom<-1

### set sqrt(2) for the "zero.removal" option ###

    if(correction.method == "zero.removal") denom<-sqrt(2)

### perform the orthogonal transformation ###

    replace1<-(sub.mean.diff[mid2]-sub.mean.diff[mid])/denom
    temp2<-sub.mean.diff[c(-mid,-mid2)]
    temp2<-c(temp2,replace1)
   }

### collect the transformed variables into the vector ###

   temp<-c(temp,temp2)
  }

  resp.mean<-abs(temp)

### multiply the correction factor for the "zero.correction" option ###

  if(correction.method == "zero.correction")
  {
   correction<-sqrt((ngroup-1)/ngroup)
   correction<-correction[-endpos]
   resp.mean<-resp.mean*correction
  }
   
### update variables used for the contrast vector calculation ###

  group<-group[-endpos]
  n<-n-1
  r.means <- tapply(resp.mean,group,mean)
 }

### set correction.method to be "none" if specified other than those in the option ###

 else
 {
  correction.method <- "none"
  n <- tapply(y, group, length)
  resp.mean <- abs(y - means[group])
  r.means <- tapply(resp.mean,group,mean)
 }

### calculate contrast vector (followed Neuhauser and Hothorn (2000)) ###

 k <- length(n)
 a1<-double(k)
 a1[1]<--1
 end<-k-1

 for(i in 2:end)
 {
  start<-k-i+1
  finish<-k-1
  a1[i]<-sum(1/c(start:finish))-1
 }

 a1[k]<-sum(1/c(1:end))
 a2<--a1[k:1]

### calculate the statistic ###

 s <- sqrt(sum((resp.mean - r.means[group])^2)/(sum(n)-k))
 T1 <- sum(r.means*a1)/(s*sqrt(sum(a1^2/n)))
 T2 <- sum(r.means*a2)/(s*sqrt(sum(a2^2/n)))
 statistic <- max(T1,T2)

### calculate correlation matrix calculation (followed Bechhofer and Dunnet (1982)) ###

 rho <- sum(a1*a2/n)/sqrt(sum(a1^2/n)*sum(a2^2/n))
 corr <- matrix(data=c(1,rho,rho,1),nrow=2,ncol=2)
 df <- sum(n)-k

### calculate the p-value ###

 p.value <- 1 - pmvt(lower=c(-Inf,-Inf),upper=c(statistic,statistic),df=df,corr=corr)[1]

 if(tail=="right") 
 {
  METHOD<-paste(METHOD, "(right-tailed)")
  p.value <- 1 - pmvt(lower=c(-Inf,-Inf),upper=c(statistic,statistic),df=df,corr=corr)[1]
 }

 else if(tail=="left") 
 {
  METHOD<-paste(METHOD, "(left-tailed)")
  p.value <- 1 - pmvt(lower=c(statistic,statistic),upper=c(Inf,Inf),df=df,corr=corr)[1]
 }

 else
 {
  tail = "both"
  METHOD<-paste(METHOD, "(two-tailed)")
  p.value <- 1 - pmvt(lower=c(-abs(statistic),-abs(statistic)),upper=c(abs(statistic),abs(statistic)),df=df,corr=corr)[1]
 }

### store the non-boostrap p-value ###

 non.bootstrap.p.value <- p.value

### perform bootstrapping (followed Lim and Loh(1996)) ###

 if(bootstrap==TRUE) 
 {
  METHOD = paste("bootstrap", METHOD)

### re-initialize some of the variables ###

   n<-original.n
   group<-original.group

### step 2 of Lim and Loh (1996): initialize variables ###

  R<-0
  N<-length(y)

### step 3 of Lim and Loh (1996): calculate the fractional trimmed mean ###

  frac.trim.alpha = 0.2
  b.trimmed.mean<-function(y)
  {
   nn<-length(y)
   wt<-rep(0,nn)
   y2<-y[order(y)]
   lower<-ceiling(nn*frac.trim.alpha)+1
   upper<-floor(nn*(1-frac.trim.alpha))

   if(lower>upper) stop("frac.trim.alpha value is too large")
   m<-upper-lower+1
   frac<-(nn*(1-2*frac.trim.alpha)-m)/2
   wt[lower-1]<-frac
   wt[upper+1]<-frac
   wt[lower:upper]<-1
   return(weighted.mean(y2,wt))
  }
  
  b.trim.means <- tapply(y, group, b.trimmed.mean)
  rm<-y-b.trim.means[group]

### step 7 of Lim and Loh (1996): enter a loop ###

  for (j in 1:num.bootstrap)
  {

### re-initialize some of the variables ###

   n<-original.n
   group<-original.group

### step 4 of Lim and Loh (1996): obtain a bootstrap sample ###

   sam<-sample(rm,replace=TRUE)
   boot.sample <- sam

### step 5 of Lim and Loh (1996): smooth the variables if n_i < 10 for at least one sample size ###

   if(min(n) < 10)
   {
    U<-runif(1)-0.5
    means <- tapply(y, group, mean)
    v <- sqrt(sum((y - means[group])^2)/N)
    boot.sample <- ((12/13)^(0.5))*(sam + v*U)
   }

### step 6 of Lim and Loh (1996): compute the bootstrap statistic, and increment R to R + 1 if necessary ###

   if(location=="mean")
   {
    boot.means <- tapply(boot.sample, group, mean)
   }

   else if(location=="median")
   { 
    boot.means <- tapply(boot.sample, group, median)
   }

   else 
   {location="trim.mean";
    trimmed.mean.2<-function(boot.sample) mean(boot.sample, trim=trim.alpha) 
    boot.means <- tapply(boot.sample, group, trimmed.mean.2)
   }

### calculate bootstrap statistic ###

   resp.boot.mean <- abs(boot.sample - boot.means[group])

### multiply the correction factor to each observation if "correction.factor" is chosen ###

   if(correction.method == "correction.factor") 
   {
    correction<-1/sqrt(1-1/ngroup)
    resp.boot.mean<-resp.boot.mean*correction
    r.boot.means <- tapply(resp.boot.mean,group,mean)
   }

### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
### or "zero.correction" (Noguchi and Gel, 2009).                            ###

   else if(correction.method == "zero.correction"||correction.method == "zero.removal")
   {

### set up variables for calculating the deviation from center ###

    mean.diff <- boot.sample-boot.means[group]
    k<-length(n)
    temp<-double()
    endpos<-double()
    startpos<-double()

### calculate the absolute deviation from mean and remove zeros ###

    for(i in 1:k)
    {
     group.size<-n[i]
     j<-i-1

### calculate the starting and ending index of each group ###

     if(i==1)start<-1

     else start<-sum(n[1:j])+1

     startpos<-c(startpos,start)
     end<-sum(n[1:i])
     endpos<-c(endpos,end)

### extract the deviation from center for the ith group ###

     sub.mean.diff<-mean.diff[start:end]
     sub.mean.diff<-sub.mean.diff[order(sub.mean.diff)]

### remove structural zero for the odd-sized group ###

     if(group.size%%2==1)
     {
      mid<-(group.size+1)/2
      temp2<-sub.mean.diff[-mid]
     }

### remove structural zero for the even-sized group ###

     if(group.size%%2==0)
     {
      mid<-group.size/2
      mid2<-mid+1

### set up the denominator value for the transformation ###

### set 1 for the "zero.correction" option ###

      denom<-1

### set sqrt(2) for the "zero.removal" option ###

      if(correction.method == "zero.removal") denom<-sqrt(2)

### perform the orthogonal transformation ###

      replace1<-(sub.mean.diff[mid2]-sub.mean.diff[mid])/denom
      temp2<-sub.mean.diff[c(-mid,-mid2)]
      temp2<-c(temp2,replace1)
     }

### collect the transformed variables into the vector ###

     temp<-c(temp,temp2)
    }

    resp.boot.mean<-abs(temp)

### multiply the correction factor for the "zero.correction" option ###

    if(correction.method == "zero.correction")
    {
     correction<-sqrt((ngroup-1)/ngroup)
     correction<-correction[-endpos]
     resp.boot.mean<-resp.boot.mean*correction
    }
   
### update variables used for the contrast vector calculation ###

    group<-group[-endpos]
    n<-n-1
    r.boot.means <- tapply(resp.boot.mean,group,mean)
   }

   else
   {
    r.boot.means <- tapply(resp.boot.mean,group,mean)
   }

   boot.s <- sqrt(sum((resp.boot.mean - r.boot.means[group])^2)/(sum(n)-k)) 
   boot.T1 <- sum(r.boot.means*a1)/(boot.s*sqrt(sum(a1^2/n)))
   boot.T2 <- sum(r.boot.means*a2)/(boot.s*sqrt(sum(a2^2/n)))
   statistic2 <- max(boot.T1,boot.T2)

   if(tail=="right")
   {
    if(statistic2 > statistic) R<-R+1
   }

   else if(tail=="left")
   {
    if(statistic2 < statistic) R<-R+1
   }

   else
   {tail="both";
    if(abs(statistic2) > abs(statistic)) R<-R+1
   }
  }
   
### step 8 of Lim and Loh (1996): calculate the bootstrap p-value ###

  p.value <- R/num.bootstrap
 }

### display output ###

 STATISTIC=statistic
 names(STATISTIC)="Test Statistic"
 
 structure(list(statistic = STATISTIC, p.value=p.value, method = METHOD,
 data.name = DNAME, non.bootstrap.p.value=non.bootstrap.p.value), class = "htest")
}

