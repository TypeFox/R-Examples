#
# vim:set ff=unix expandtab ts=2 sw=2:
plotC14Pool<-structure(
  function #Plots the output of \code{\link{getC14}} for each pool over time
  ### This function produces a plot with the Delta14C in the atmosphere and the Delta14C of each pool obtained after a call to \code{\link{getC14}}.
    (t,  ##<< A vector containing the time points for plotting.
     mat, ##<< A matrix object obtained after a call to \code{\link{getC14}}
     inputFc, ##<< A Data Frame object containing values of atmospheric Delta14C per time. First column must be time values, second column must be Delta14C values in per mil.
     col,   ##<< A color palette specifying color lines for each pool (columns of \code{mat}).
     ...    ##<< Other arguments passed to \code{plot}.
     )
    {
     n=dim(mat)[2]
     plot(inputFc,type="l",...)
         for(i in 1:n){
         lines(t,mat[,i],col=col[i])
         }
    }
  ,
  ex=function(){
    years=seq(1901,2009,by=0.5)
    LitterInput=700 
    
    Ex=ThreepFeedbackModel14(
      t=years,ks=c(k1=1/2.8, k2=1/35, k3=1/100),
      C0=c(200,5000,500), F0=c(0,0,0),
      In=LitterInput, a21=0.1,a12=0.01,a32=0.005,a23=0.001,inputFc=C14Atm_NH
    )
    C14t=getF14(Ex)
    
    pal=rainbow(3)
    plotC14Pool(
      t=years,mat=C14t,inputFc=C14Atm_NH,
      col=pal,xlab="Time (yrs)",
      ylab="Delta 14C (per mil)",
      xlim=c(1950,2000)
    )
  }
)
