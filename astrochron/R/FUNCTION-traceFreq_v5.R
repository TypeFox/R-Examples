### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### traceFreq: trace frequency drift using graphical interface 
###                      (SRM: June 12-14, 2013; June 16, 2013; June 26, 2013
###                            December 8, 2013; January 16, 2015)
###  
###########################################################################

traceFreq <- function (spec,color=2,h=6,w=4,ydir=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,ncolors=100,pl=1,ln=F)
{

    cat("\n----            FREQUENCY DOMAIN MINIMAL TUNING:              ----\n")
    cat("---- INTERACTIVELY TRACE SPATIAL FREQUENCY DRIFT IN EHA PLOTS ----\n")
    cat("\n           *****  Select path by clicking     *****\n")
    cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")

# use fields library for access to 'tim.colors'
  spec=data.frame(spec)

# assign frequencies from first column of spec
  freq=spec[,1]

  cols=length(spec)
# assign locations for each spectrum (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(spec[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(spec[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(spec[(fixloc+1)]),start=3,stop=100))}
# assign specta
  sp=as.matrix( spec[2:cols] )

  if(ln) sp=log(sp)
  if(is.null(xmin)) xmin = min(freq)
  if(is.null(xmax)) xmax = max(freq)
  if(is.null(ymin)) ymin = min(loc)
  if(is.null(ymax)) ymax = max(loc)
  if(pl==1) pltype="o"
  if(pl==2) pltype="l"
  if(pl==3) pltype="p"

# set up plot
  dev.new(height=h,width=w)
  par(mfrow=c(1,1))
  xlimset=c(xmin,xmax)

  if (ydir == -1) 
    {
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
       ylimset=c(ymax,ymin)
       image(freq,loc,sp,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),xlab="Frequency",ylab="Location",main="Frequency Domain Minimal Tuning")
     }

    if (ydir == 1) 
     {
# useRaster=T results in a faster plotting time.
       ylimset=c(ymin,ymax)
       image(freq,loc,sp,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),useRaster=T,xlab="Frequency",ylab="Location",main="Frequency Domain Minimal Tuning")       
   }

# Now overlay x-y plot for graphical interface
       par(new=T)
       plot(-1,-1,xlim=xlimset,ylim=ylimset,xaxs="i",yaxs="i",yaxt='n',bty='n',ylab="",xlab="")
# transparent black
       if (color == 1) setcolor="#00000070"
# transparent white: alpha (transparency) of 210 in hexadecimal = D2       
       if (color == 2) setcolor="#FFFFFFD2"
# transparent yellow
       if (color == 3) setcolor="#FFFF00D2"
       ff = locator(n = length(loc), type = pltype ,col=setcolor,lwd=3)

       out = data.frame(cbind(ff$y,ff$x))
       colnames(out)[1] = 'Depth/Height'
       colnames(out)[2] = 'Frequency'
       
# Sort to ensure Depth/Height is in increasing order
       out <- out[order(out[1],na.last=NA,decreasing=F),]

       cat("\n * Tracing complete.\n")
       cat("\n * Now:\n")
       cat("   (1) Convert the spatial frequencies to sedimentation rates using function 'freq2sedrate'\n")
       cat("   (2) Integrate sedimentation rate curve using function 'sedrate2time'\n")
       cat("   (3) Tune proxy record using function 'tune'\n")

       return(out)

### END function traceFreq
}
