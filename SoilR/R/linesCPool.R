#
# vim:set ff=unix expandtab ts=2 sw=2:
linesCPool<-structure(
  function #Add lines with the output of \code{\link{getC14}}, \code{\link{getC}}, or \code{\link{getReleaseFlux}} to an existing plot
  ### This function adds lines to a plot with the C content, the C release, or Delta 14C value of each pool over time. Needs as input a matrix obtained after a call to \code{\link{getC14}}, \code{\link{getC}}, or \code{\link{getReleaseFlux}}.
  (t,  ##<< A vector containing the time points for plotting.
   mat, ##<< A matrix object obtained after a call to \code{\link{getC14}}, \code{\link{getC}}, or \code{\link{getReleaseFlux}}.
   col,   ##<< A color palette specifying color lines for each pool (columns of \code{mat}).
   ...    ##<< Other arguments passed to \code{plot}.
  )
  {
      n=dim(mat)[2]
      for(i in 1:n){
        lines(t,mat[,i],col=col[i],...)
      }
  }
  ,
  ex=function(){
    years=seq(1901,2009,by=0.5)
    LitterInput=700 
    
    Ex=ThreepFeedbackModel14(t=years,ks=c(k1=1/2.8, k2=1/35, k3=1/100),
                             C0=c(200,5000,500),F0=c(0,0,0), In=LitterInput, 
                             a21=0.1,a12=0.01,a32=0.005,a23=0.001,inputFc=C14Atm_NH)
    Ct=getC(Ex)
    
    pal=rainbow(3)
    plotCPool(t=years,mat=Ct,col=pal,xlab="Time (yrs)",
              ylab="Carbon stocks",ylim=c(min(Ct),max(Ct)))
    
    LitterInput2=350 
    
    Ex2=ThreepFeedbackModel14(t=years,ks=c(k1=1/2.8, k2=1/35, k3=1/100),
                              C0=c(200,5000,500),F0=c(0,0,0), In=LitterInput2, 
                              a21=0.1,a12=0.01,a32=0.005,a23=0.001,inputFc=C14Atm_NH)
    Ct2=getC(Ex2)
    
    linesCPool(t=years,mat=Ct2,col=pal,lwd=2)
    
  }
)
