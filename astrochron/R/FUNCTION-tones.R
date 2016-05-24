### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### tones: determine all possible difference and combinations tones
###        from a set of frequencies, find the closest one to a
###        specified frequency (SRM: January 19, 2014; January 24, 2014)
###
###########################################################################

tones <- function (a=NULL,freqs=NULL,f=T)
{
 
cat("\n----- COMPUTING DIFFERENCE AND COMBINATION TONES -----\n")

# if not defined, set freqs to default
# a=etp(dt=5,esinw=T,tmax=2000)
# eha(a,win=100000,pad=20000,sigID=T)
if(is.null(freqs)) 
 {
    cat("\n**** WARNING: using default astronomical frequencies from 0-5 Ma, LA04.\n")
    freqs=c(1/404.8583,1/126.2626,1/96.06148,1/53.73455,1/41.13534,1/39.3391,1/28.82675,1/23.67985,1/22.37136,1/18.99696)
 }

freqs=data.frame(freqs)
nfreq=length(freqs[,1])

# sort freqs into increasing order
freqs=data.frame(sort(freqs[,1]))

if(f) cat("\n ** Difference tones (cycles/ka):")
if(!f) cat("\n ** Difference tones (ka):")

tone <- rep(NA,nfreq*nfreq)
tone <- as.double(tone)
toneType <- character(nfreq*nfreq)
k=1
for (i in 1:nfreq)
  {
    for (j in i:nfreq)
      {
        if(i != j)
          {
           tone[k] = freqs[j,1]-freqs[i,1]
           toneType[k] = "-"
           if(f) cat("\n",k,".",freqs[j,1],"-",freqs[i,1],"=",tone[k])
           if(!f) cat("\n",k,".",1/freqs[j,1],"-",1/freqs[i,1],"=",1/tone[k])
           k = k + 1
          }
      }
  }

cat("\n")

if(f) cat("\n ** Combination tones (cycles/ka):")
if(!f) cat("\n ** Combination tones (ka):")

for (i in 1:nfreq)
  {
    for (j in i:nfreq)
      {
         if(i != j) 
           {
             tone[k] = freqs[j,1]+freqs[i,1]
             toneType[k] = "+"
             if(f) cat("\n",k,".",freqs[j,1],"+",freqs[i,1],"=",tone[k])
             if(!f) cat("\n",k,".",1/freqs[j,1],"+",1/freqs[i,1],"=",1/tone[k])
             k = k + 1
            } 
      }
  }

# sort into increasing order, remove NAs, report index
ii=sort.int(tone,na.last=NA,index.return=T)
out=tone[ii$ix]
iout=length(out)

par(mfrow=c(1,1))
plot(out,main="All Combination and Difference Tones",ylab="Frequency (cycles/ka)",cex=1.2)
text(1:iout,out,labels=ii$ix,cex=0.35,col="red")
if(!is.null(a)) 
 {
   abline(h=a,col="red",lty=2)
   text(1,a,labels="Your frequency",col="red",adj=0)
 }

if(!is.null(a))
   {
# note: this tone still includes NAs   
# does not check to see if multiple tones yield best fit!
     minout=abs(a-tone[1])
     minid=1
     for (i in 2:iout)
       {
         test=abs(a-tone[i])
          if(test < minout)
            {
               minout=test
               minid=i
             }  
        }     
      cat("\n \n Best fit identified for tone number",minid)   
    }  

### END function tones
}
