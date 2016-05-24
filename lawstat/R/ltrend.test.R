`ltrend.test` <-
function (y, group, score=NULL, location = c("median", "mean", "trim.mean"), 
tail = c("right","left","both"), trim.alpha = 0.25, bootstrap = FALSE, num.bootstrap = 1000, 
correction.method = c("none","correction.factor","zero.removal","zero.correction"), 
correlation.method = c("pearson","kendall","spearman"))
{

### assign score to each group ###

 if (is.null(score)) 
 {
  score<-group
 }

### stop the code if the length of y does not match the length of group ###

 if (length(y)!=length(group))
 {
  stop("the length of the data (y) does not match the length of the group")
 }

### assign location, tail, and a correction method ###

 location<-match.arg(location)
 tail<-match.arg(tail)
 correlation.method<-match.arg(correlation.method)
 correction.method<-match.arg(correction.method)
 DNAME<-deparse(substitute(y))
 y<-y[!is.na(y)]
 score<-score[!is.na(y)]
 group<-group[!is.na(y)]
 
### stop the code if the location "trim.mean" is selected and trim.alpha is too large ###

 if ((location == "trim.mean") & (trim.alpha == 1)) 
 {
  stop("trim.alpha value of 0 to 0.5 should be provided for the trim.mean location")
 }
 
### sort the order just in case the input is not sorted by group ###
 
 reorder<-order(group)
 group<-group[reorder]
 y<-y[reorder]
 score<-score[reorder]
 gr<-score
 group<-as.factor(group)

### define the measure of central tendency (mean, median, trimmed mean) ###
 
 if (location == "mean") 
 {
  means<-tapply(y, group, mean)
  METHOD<-"ltrend test based on classical Levene's procedure using the group means"
 }
 
 else if (location == "median") 
 {
  means<-tapply(y, group, median)
  METHOD<-"ltrend test based on the modified Brown-Forsythe Levene-type procedure using the group medians"
 }
 
 else 
 {
  location<-"trim.mean"
  trimmed.mean<-function(y) mean(y, trim = trim.alpha)
  means<-tapply(y, group, trimmed.mean)
  METHOD<-"ltrend test based on the modified Brown-Forsythe Levene-type procedure using the group trimmed means"
 }

### calculate the sample size of each group and absolute deviation from center ###   
 
 n<-tapply(y, group, length)
 ngroup<-n[group]
 resp.mean<-abs(y - means[group])

### assign no correction technique if the central tendency is median, and ###
### any technique other than "correction.factor" is chosen                ###
 
 if(location != "median" && correction.method != "correction.factor")
 {
  METHOD<-paste(METHOD,"(",correction.method,"not applied because the location is not set to median",")")
  correction.method<-"none"
 }
 
### multiply the correction factor to each observation if "correction.factor" is chosen ###

 if(correction.method == "correction.factor") 
 {
  METHOD<-paste(METHOD,"with correction factor applied")
  correction<-1/sqrt(1-1/ngroup)
  resp.mean<-resp.mean*correction
 }
 
### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
### or "zero.correction" (Noguchi and Gel, 2009).                            ###
 
 if(correction.method =="zero.removal"||correction.method=="zero.correction")
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
  
  resp.mean <- y - means[group]
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

   sub.resp.mean<-resp.mean[start:end]
   sub.resp.mean<-sub.resp.mean[order(sub.resp.mean)]

### remove structural zero for the odd-sized group ###
   
   if(group.size%%2==1)
   {
    mid<-(group.size+1)/2
    temp2<-sub.resp.mean[-mid]
    
### multiply the correction factor for the "zero.correction" option ###

    if(correction.method=="zero.correction")
    {
     ntemp<-length(temp2)+1
     correction<-sqrt((ntemp-1)/ntemp)
     temp2<-correction*temp2
    }
   }

### remove structural zero for the even-sized group ###
   
   if(group.size%%2==0)
   {
    mid<-group.size/2
    
### set up the denominator value for the transformation ###

### set sqrt(2) for the "zero.removal" option ###

    if(correction.method=="zero.removal") 
    {
     denom<-sqrt(2)
    }

### set 1 for the "zero.correction" option ###
    
    else
    {
     denom<-1
    }
    
### perform the orthogonal transformation ###

    replace1<-(sub.resp.mean[mid+1]-sub.resp.mean[mid])/denom
    temp2<-sub.resp.mean[c(-mid,-mid-1)]
    temp2<-c(temp2,replace1)
    
### multiply the correction factor for the "zero.correction" option ###

    if(correction.method=="zero.correction")
    {
     ntemp<-length(temp2)+1
     correction<-sqrt((ntemp-1)/ntemp)
     temp2<-correction*temp2
    }
   }

### collect the transformed variables into the vector ###

   temp<-c(temp,temp2)
  }
  
### calculate the absolute deviation from center with modifications ###

  ngroup2<-ngroup[-endpos]-1
  resp.mean<-abs(temp)
  zero.removal.gr<-gr[-endpos]
 }
 
### set correction.method to be "none" if specified other than those in the option ###

 else
 {
  correction.method="none"
 }
 
### calculate z ###

 mu<-mean(resp.mean)
 z<-as.vector(resp.mean - mu)
 
### set zero.removal.gr as d for methods with structural zero removal ###

 if(correction.method=="zero.removal"||correction.method=="zero.correction")
 {
  d<-as.numeric(zero.removal.gr)
 }
 
### set the original gr as d otherwise ###

 else
 {
  d<-as.numeric(gr)
 }
 
### calculate the t-statistic using a simple linear regression ###

 t.statistic<-summary(lm(z ~ d))$coefficients[2, 3]
 df<-summary(lm(z ~ d))$df[2]
 
### calculate if the tail is set to "left" ###

 if (correlation.method == "pearson")
 {
### calculate the correlation between d and z ###

  correlation<-cor(d,z,method="pearson")

  if (tail == "left") 
  {
   METHOD<-paste(METHOD, "(left-tailed with Pearson correlation coefficient)")
   p.value<-pt(t.statistic,df,lower.tail=TRUE)
   log.p.value<-pt(t.statistic,df,lower.tail=TRUE,log.p=TRUE)
   log.q.value<-pt(t.statistic,df,lower.tail=FALSE,log.p=TRUE)
  }
 
### calculate if the tail is set to "right" ###

  else if (tail == "right") 
  {
   METHOD<-paste(METHOD, "(right-tailed with Pearson correlation coefficient)")
   p.value<-pt(t.statistic,df,lower.tail=FALSE)
   log.p.value<-pt(t.statistic,df,lower.tail=FALSE,log.p=TRUE)
   log.q.value<-pt(t.statistic,df,lower.tail=TRUE,log.p=TRUE)
  }
 
### calculate if the tail is set to "both" ###

  else 
  {tail<-"both";
   METHOD<-paste(METHOD, "(two-tailed with Pearson correlation coefficient)")
   p.value<-pt(t.statistic,df,lower.tail=TRUE)
   log.p.value<-pt(t.statistic,df,lower.tail=TRUE,log.p=TRUE)
   log.q.value<-pt(t.statistic,df,lower.tail=FALSE,log.p=TRUE)
  
   if(p.value >= 0.5) 
   {
    p.value<-pt(t.statistic,df,lower.tail=FALSE)
    log.p.value<-pt(t.statistic,df,lower.tail=FALSE,log.p=TRUE)
    log.q.value<-pt(t.statistic,df,lower.tail=TRUE,log.p=TRUE)
   }
  
   p.value<-p.value*2
   log.p.value<-log.p.value+log(2)
   log.q.value<-log(1-p.value)
  }
 }
 
 else if (correlation.method == "kendall")
 {
### calculate the correlation between d and z ###

  correlation<-cor(d,z,method="kendall")

  if (tail == "left") 
  {
   METHOD<-paste(METHOD, "(left-tailed with Kendall correlation coefficient)")
   p.value.temp<-Kendall(d,z)$sl

   if(correlation < 0)
   {
    p.value<-p.value.temp/2
   }

   else
   {
    p.value<-1-p.value.temp/2
   }
   q.value<-1-p.value
   log.p.value<-log(p.value)
   log.q.value<-log(q.value)
  }

  if (tail == "right") 
  {
   METHOD<-paste(METHOD, "(right-tailed with Kendall correlation coefficient)")
   p.value.temp<-Kendall(d,z)$sl

   if(correlation > 0)
   {
    p.value<-p.value.temp/2
   }

   else
   {
    p.value<-1-p.value.temp/2
   }
   q.value<-1-p.value  
   log.p.value<-log(p.value)
   log.q.value<-log(q.value)
  }

  if (tail == "both") 
  {
   METHOD<-paste(METHOD, "(two-tailed with Kendall correlation coefficient)")
   p.value<-Kendall(d,z)$sl
   q.value<-1-p.value 
   log.p.value<-log(p.value)
   log.q.value<-log(q.value)
  }
 }

 else
 {
### calculate the correlation between d and z ###

  correlation<-cor(d,z,method="spearman")

  if (tail == "left") 
  {
   METHOD<-paste(METHOD, "(left-tailed with Spearman correlation coefficient)")
   p.value.temp<-spearman.test(d,z)[5]
   if(correlation < 0)
   {
    p.value<-p.value.temp/2
   }

   else
   {
    p.value<-1-p.value.temp/2
   }
   
   q.value<-1-p.value   
   log.p.value<-log(p.value)
   log.q.value<-log(q.value)
  }

  if (tail == "right") 
  {
   METHOD<-paste(METHOD, "(right-tailed with Spearman correlation coefficient)")
   p.value.temp<-spearman.test(d,z)[5]
   if(correlation > 0)
   {
    p.value<-p.value.temp/2
   }

   else
   {
    p.value<-1-p.value.temp/2
   }
   
   q.value<-1-p.value 
   log.p.value<-log(p.value)
   log.q.value<-log(q.value)
  }

  if (tail == "both") 
  {
   METHOD<-paste(METHOD, "(two-tailed with Spearman correlation coefficient)")
   p.value<-spearman.test(d,z)[5]
   q.value<-1-p.value 
   log.p.value<-log(p.value)
   log.q.value<-log(q.value)
  }
 }
 
### store the non-boostrap p-value ###

 non.bootstrap.p.value <- p.value
 
### perform bootstrapping (followed Lim and Loh(1996)) ###

 if(bootstrap == TRUE)
 {
  METHOD = paste("bootstrap", METHOD)
  
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
  
  b.trim.means<-tapply(y, group, b.trimmed.mean)
  rm<-y-b.trim.means[group]
  
### step 7 of Lim and Loh (1996): enter a loop ###

  for (j in 1:num.bootstrap)
  {

### step 4 of Lim and Loh (1996): obtain a bootstrap sample ###

   sam<-sample(rm,replace=TRUE)
   boot.sample<-sam
   
### step 5 of Lim and Loh (1996): smooth the variables if n_i < 10 for at least one sample size ###

   if(min(n)<10)
   {
    U<-runif(1)-0.5
    means<-tapply(y, group, mean)
    v<-sqrt(sum((y-means[group])^2)/N)
    boot.sample<-((12/13)^(0.5))*(sam + v*U)
   }
   
### step 6 of Lim and Loh (1996): compute the bootstrap statistic, and increment R to R + 1 if necessary ###
   
   if(location=="mean")
   {
    boot.means<-tapply(boot.sample, group, mean)
   }
   
   else if(location=="median")
   {
    boot.means<-tapply(boot.sample, group, median)
   }
   
   else 
   {location<-"trim.mean";
    trimmed.mean.2<-function(boot.sample) mean(boot.sample, trim=trim.alpha) 
    boot.means<-tapply(boot.sample, group, trimmed.mean.2)
   }

### calculate bootstrap statistic ###

   resp.boot.mean<-abs(boot.sample - boot.means[group])

   if(correction.method == "correction.factor") 
   {
    correction<-1/sqrt(1-1/ngroup)
    resp.boot.mean<-resp.boot.mean*correction
   }

### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
### or "zero.correction" (Noguchi and Gel, 2009).                            ###
 
   if(correction.method =="zero.removal"||correction.method=="zero.correction")
   {

### set up variables for calculating the deviation from center ###
  
    resp.mean <- boot.sample - boot.means[group]
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

     sub.resp.mean<-resp.mean[start:end]
     sub.resp.mean<-sub.resp.mean[order(sub.resp.mean)]

### remove structural zero for the odd-sized group ###
   
     if(group.size%%2==1)
     {
      mid<-(group.size+1)/2
      temp2<-sub.resp.mean[-mid]
    
### multiply the correction factor for the "zero.correction" option ###

      if(correction.method=="zero.correction")
      {
       ntemp<-length(temp2)+1
       correction<-sqrt((ntemp-1)/ntemp)
       temp2<-correction*temp2
      }
     }

### remove structural zero for the even-sized group ###
   
      if(group.size%%2==0)
      {
       mid<-group.size/2
    
### set up the denominator value for the transformation ###

### set sqrt(2) for the "zero.removal" option ###

      if(correction.method=="zero.removal") 
      {
       denom<-sqrt(2)
      }

### set 1 for the "zero.correction" option ###
    
      else
      {
       denom<-1
      }
    
### perform the orthogonal transformation ###

      replace1<-(sub.resp.mean[mid+1]-sub.resp.mean[mid])/denom
      temp2<-sub.resp.mean[c(-mid,-mid-1)]
      temp2<-c(temp2,replace1)
    
### multiply the correction factor for the "zero.correction" option ###

      if(correction.method=="zero.correction")
      {
       ntemp<-length(temp2)+1
       correction<-sqrt((ntemp-1)/ntemp)
       temp2<-correction*temp2
      }
     }

### collect the transformed variables into the vector ###

     temp<-c(temp,temp2)
    }
  
### calculate the absolute deviation from center with modifications ###

    ngroup2<-ngroup[-endpos]-1
    resp.boot.mean<-abs(temp)
    zero.removal.gr<-gr[-endpos]
   }

### set zero.removal.gr as d for methods with structural zero removal ###

   if(correction.method=="zero.removal"||correction.method=="zero.correction")
   {
    d<-as.numeric(zero.removal.gr)
   }
 
### set the original gr as d otherwise ###

   else
   {
    d<-as.numeric(gr)
   }

   boot.mu<-mean(resp.boot.mean)
   boot.z<-as.vector(resp.boot.mean-boot.mu)
   correlation2<-cor(boot.z,d,method=correlation.method)
   
   if(tail=="right")
   {
    if(correlation2>correlation) R<-R+1
   }
   
   else if(tail=="left")
   {
    if(correlation2<correlation) R<-R+1
   }
   
   else
   {tail="both";
    if(abs(correlation2)>abs(correlation)) R<-R+1
   }
  }
  
### step 8 of Lim and Loh (1996): calculate the bootstrap p-value ###

  p.value <- R/num.bootstrap
 }
 
### display output ###

 STATISTIC = correlation
 
 names(STATISTIC) = "Test Statistic (Correlation)"
 structure(list(statistic = STATISTIC, p.value = p.value, method = METHOD, 
 data.name = DNAME, t.statistic = t.statistic, non.bootstrap.p.value = non.bootstrap.p.value, 
 log.p.value = log.p.value, log.q.value = log.q.value), class = "htest")
}

