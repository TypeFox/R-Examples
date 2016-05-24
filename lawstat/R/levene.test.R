`levene.test` <-
function(y, group, location=c("median", "mean", "trim.mean"), trim.alpha=0.25,
bootstrap = FALSE, num.bootstrap=1000, kruskal.test=FALSE, 
correction.method=c("none","correction.factor","zero.removal","zero.correction"))
{
 
### stop the code if the length of y does not match the length of group ###

 if (length(y)!=length(group))
 {
  stop("the length of the data (y) does not match the length of the group")
 }

### assign location, tail, and a correction method ###

 location<-match.arg(location)
 correction.method<-match.arg(correction.method)
 DNAME = deparse(substitute(y))
 y<-y[!is.na(y)]
 group<-group[!is.na(y)]

### stop the code if the location "trim.mean" is selected and trim.alpha is too large ###

 if ((location=="trim.mean")&(trim.alpha==1))
 {
  stop("trim.alpha value of 0 to 0.5 should be provided for the trim.mean location")
 }
 
### sort the order just in case the input is not sorted by group ###
 
 reorder<-order(group)
 group<-group[reorder]
 y<-y[reorder]
 gr<-group
 group<-as.factor(group) # precautionary
 

### define the measure of central tendency (mean, median, trimmed mean) ###
 
 if(location=="mean")
 {
  means<-tapply(y, group, mean)
  METHOD<-"classical Levene's test based on the absolute deviations from the mean"
 }
 
 else if(location=="median")
 {
  means<-tapply(y, group, median)
  METHOD="modified robust Brown-Forsythe Levene-type test based on the absolute deviations from the median"
 }
 
 else 
 {location="trim.mean";
  trimmed.mean<-function(y) mean(y, trim=trim.alpha)
  means<-tapply(y, group, trimmed.mean)
  METHOD<-"modified robust Levene-type test based on the absolute deviations from the trimmed mean"
 }

### calculate the sample size of each group and absolute deviation from center ### 
 
 n<-tapply(y, group, length)
 resp.mean<-abs(y - means[group])
 ngroup<-n[group]

### assign no correction technique if the central tendency is median, and ###
### any technique other than "correction.factor" is chosen                ###
 
 if(location!="median" && correction.method!="correction.factor")
 {
  METHOD<-paste(METHOD,"(",correction.method,"not applied because the location is not set to median",")")
  correction.method<-"none"
 }
 
### multiply the correction factor to each observation if "correction.factor" is chosen ###

 if(correction.method=="correction.factor")
 {
  METHOD<-paste(METHOD,"with correction factor")
  correction<-sqrt(ngroup/(ngroup-1))
  resp.mean<-correction*resp.mean
 }
  
### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
### or "zero.correction" (Noguchi and Gel, 2009).                            ###

 if(correction.method=="zero.removal"||correction.method=="zero.correction")
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

  resp.mean<-y-means[group]
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

  resp.mean<-abs(temp)
  zero.removal.group<-group[-endpos]
 }
 
### set correction.method to be "none" if specified other than those in the option ###

 else
 {
  correction.method="none"
 }
 
### calculate statistic and p-value when structural zeros are removed ###

### set d depending on whether structural zero removal method is used ###
 
 if(correction.method=="zero.removal"||correction.method=="zero.correction")
 {
  d<-zero.removal.group
 }
 else
 {
  d<-group
 }

### if the Kruskal-Wallis test is not used ###

 if(kruskal.test == FALSE)
 {
  statistic<-anova(lm(resp.mean ~ d))[1, 4]
  p.value<-anova(lm(resp.mean ~ d))[1, 5]
 }

### if the Kruskal-Wallis test is used ###

 else
 {
  METHOD<-paste("rank-based (Kruskal-Wallis)", METHOD)
  ktest<-kruskal.test(resp.mean,d)
  statistic<-ktest$statistic
  p.value=ktest$p.value
 }
  
### store the non-boostrap p-value ###

 non.bootstrap.p.value<-p.value
 
### perform bootstrapping (followed Lim and Loh(1996)) ###

 if(bootstrap == TRUE)
 {
  METHOD<-paste("bootstrap", METHOD)
  
### step 2 of Lim and Loh (1996): initialize variables ###

  R<-0
  N<-length(y)
  
### step 3 of Lim and Loh (1996): calculate the fractional trimmed mean ###

  frac.trim.alpha<-0.2
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
   boot.sample <- sam
   
### step 5 of Lim and Loh (1996): smooth the variables if n_i < 10 for at least one sample size ###

   if(min(n) < 10)
   {
    U<-runif(1)-0.5
    means<-tapply(y, group, mean)
    v<-sqrt(sum((y - means[group])^2)/N)
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
   {location="trim.mean";
    trimmed.mean.2<-function(boot.sample) mean(boot.sample, trim=trim.alpha) 
    boot.means<-tapply(boot.sample, group, trimmed.mean.2)
   }
   
### calculate bootstrap statistic ###

   resp.boot.mean<-abs(boot.sample - boot.means[group])

### multiply the correction factor to each observation if "correction.factor" is chosen ###

   if(correction.method=="correction.factor")
   {
    correction<-sqrt(ngroup/(ngroup-1))
    resp.mean<-correction*resp.boot.mean
   }

### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
### or "zero.correction" (Noguchi and Gel, 2009).                            ###

   if(correction.method=="zero.removal"||correction.method=="zero.correction")
   {

### set up variables for calculating the deviation from center ###

    resp.mean<-boot.sample-boot.means[group]
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

    resp.boot.mean<-abs(temp)
    zero.removal.group<-group[-endpos]
   }

### set d depending on whether structural zero removal method is used ###
 
   if(correction.method=="zero.removal"||correction.method=="zero.correction")
   {
    d<-zero.removal.group
   }
   else
   {
    d<-group
   }

### if the Kruskal-Wallis test is not used ###
   
   if(kruskal.test==FALSE)
   {
    statistic2 = anova(lm(resp.boot.mean ~ d))[1, 4]
   }

### if the Kruskal-Wallis test is used ###

   else
   {
   bktest<-kruskal.test(resp.boot.mean,d)
   statistic2<-bktest$statistic
   }
   
   if(statistic2 > statistic) R<-R+1
  }
  
### step 8 of Lim and Loh (1996): calculate the bootstrap p-value ###

  p.value<-R/num.bootstrap
 }
 
### display output ###

 STATISTIC=statistic
 names(STATISTIC)="Test Statistic"
 
 structure(list(statistic = STATISTIC, p.value = p.value, method = METHOD,
 data.name = DNAME, non.bootstrap.p.value = non.bootstrap.p.value), class = "htest")
}

