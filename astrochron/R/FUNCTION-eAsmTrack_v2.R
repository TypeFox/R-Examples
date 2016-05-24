### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### eAsmTrack function - (SRM: January 13, 2014; January 16, 2015; January 21, 2015)
###
### track ASM Null Hypothesis significance level minima
###########################################################################

eAsmTrack <- function (res,threshold=.5,ydir=-1,genplot=T,verbose=T)
{

 if(verbose) cat("\n----- IDENTIFYING MINIMUM Ho-SIGNIFICANCE LEVEL IN EACH WINDOW -----\n")

# ensure we have a data frame
  res=data.frame(res)

# assign sedrates from first column of res
  sedrates=res[,1]

  cols=length(res)
# assign locations for each spectrum (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(res[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(res[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(res[(fixloc+1)]),start=3,stop=100))}

# assign Ho-SL results
  Ho=as.matrix( res[2:length(res)] )

  numrec=length(loc)
  numsed=length(sedrates)
  if(verbose) cat("\n * Number of windows to analyze =",numrec,"\n")
  if(verbose) cat(" * Number of sedimentation rates =",numsed,"\n")

# loop over all windows 
# note: these are dimensioned to be much larger than needed
   HoMin<-double(numsed*numrec)
   sedMin<-double(numsed*numrec)
   locHold<-double(numsed*numrec)
   j=1
  for (i in 1:numrec)
    {
     if(verbose) cat(" * PROCESSING WINDOW=",i,"; Location=",loc[i],"\n")
# find minimum Ho-SL
     amin <- min(Ho[,i])
# identify all sedimentation rates with this Ho-SL value. This allows for the case
#  when there are multiple equivalent minima (unlike 'min' or 'which.min')
     imin=which(Ho[,i] == amin,arr.ind=TRUE)
     if(length(imin) > 1 && verbose) cat( "**** WARNING:", length(imin)," sedimentation rates have global minimum value\n")
     for (k in 1:length(imin))
       {
         HoMin[j]=Ho[imin[k],i]
         sedMin[j]=sedrates[imin[k]]
         locHold[j]=loc[i]
         j=j+1
       }
# end numrec loop
    }

   out <- data.frame(cbind(locHold[1:(j-1)],sedMin[1:(j-1)],HoMin[1:(j-1)]))
   colnames(out)<-c("Location","Sedrate","Ho-SL")
   
   out <- subset(out,(out[,3] <= threshold))

  if(genplot)
   {
     par(mfrow=c(1,2))
     if(ydir == 1) ylim=c(min(out[,1]),max(out[,1]))
     if(ydir == -1) ylim=c(max(out[,1]),min(out[,1]))
     plot(out[,2],out[,1],type="b",cex=0.75,ylim=ylim,xlab="Sedimentation Rate",ylab="Location")
     plot(out[,3],out[,1],type="b",cex=0.75,ylim=ylim,xlab="Ho-SL",ylab="Location")
   }
 
   return(out)

#### END function eAsmTrack
}
