#
# vim:set ff=unix expandtab ts=2 sw=2:
GaudinskiModel14<-structure(
  function #Implementation of a the six-pool C14 model proposed by Gaudinski et al. 2000
  ### This function creates a model as described in Gaudinski et al. 2000. 
  ### It is a wrapper for the more general functions  \code{\link{GeneralModel_14}} that can handle an arbitrary number of pools.
  ##references<< Gaudinski JB, Trumbore SE, Davidson EA, Zheng S (2000) Soil carbon cycling in a temperate forest: radiocarbon-based estimates of residence times, sequestration rates and partitioning fluxes. Biogeochemistry 51: 33-69
  (t,      ##<< A vector containing the points in time where the solution is sought. It must be specified within the same period for which the Delta 14 C of the atmosphere is provided. The default period in the provided dataset \code{\link{C14Atm_NH}} is 1900-2010.
   ks=c(kr=1/1.5,koi=1/1.5,koeal=1/4,koeah=1/80,kA1=1/3,kA2=1/75,kM=1/110),	##<< A vector of length 7 containing the decomposition rates for the 6 soil pools plus the fine-root pool. 
   C0=c(FR0=390, C10=220, C20=390, C30=1370, C40=90, C50=1800, C60=560),	##<< A vector of length 7 containing the initial amount of carbon for the 6 pools plus the fine-root pool.
   F0_Delta14C=rep(0,7), ##<< A vector of length 7 containing the initial amount of the radiocarbon fraction for the 7 pools as Delta14C values in per mil.
   LI=150,     ##<< A scalar or a data.frame object specifying the amount of litter inputs by time.
   RI=255,     ##<< A scalar or a data.frame object specifying the amount of root inputs by time.
   xi=1,   ##<< A scalar or a data.frame specifying the external (environmental and/or edaphic) effects on decomposition rates. 
   inputFc,##<< A Data Frame object containing values of atmospheric Delta14C per time. First column must be time values, second column must be Delta14C values in per mil.
   lambda=-0.0001209681, ##<< Radioactive decay constant. By default lambda=-0.0001209681 y^-1 . This has the side effect that all your time related data are treated as if the time unit was year.
   lag=0, ##<< A positive integer representing a time lag for radiocarbon to enter the system. 
   solver=deSolve.lsoda.wrapper, ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or  any other user provided function with the same interface.
 pass=FALSE  ##<< if TRUE Forces the constructor to create the model even if it is invalid 
   )	
{ 
    t_start=min(t)
    t_stop=max(t)
    if(length(ks)!=7) stop("ks must be of length = 7")
    if(length(C0)!=7) stop("the vector with initial conditions must be of length = 7")
    
    if(length(LI)==1) inputFluxes=BoundInFlux(
                                      function(t){
                                        matrix(
                                          nrow=7,ncol=1, c(RI,LI,0,0,0,0,0))
                                      },
                                      t_start,
                                      t_stop
                                      )
    if(class(LI)=="data.frame"){
      x1=LI[,1]  
      y1=LI[,2]  
      x2=RI[,1]  
      y2=RI[,2]  
      LitterFlux=function(t0){as.numeric(spline(x1,y1,xout=t0)[2])}
      RootFlux=function(t0){as.numeric(spline(x2,y2,xout=t0)[2])}
      inputFluxes= BoundInFlux(map=function(t){matrix(nrow=7,ncol=1,c(RootFlux(t),LitterFlux(t),0,0,0,0,0))}, t_start, t_stop )   
    }
    
    if(length(xi)==1) fX=function(t){xi}
    if(class(xi)=="data.frame"){
      X=xi[,1]
      Y=xi[,2]
      fX=function(t){as.numeric(spline(X,Y,xout=t)[2])}
    }
    
    A=-abs(diag(ks))
    A[3,2]=ks[2]*(98/(3+98+51))
    A[4,3]=ks[3]*(4/(94+4))
    A[6,5]=ks[5]*(24/(6+24))
    A[7,6]=ks[6]*(3/(22+3))
    A[7,2]=ks[2]*(3/(3+98+51))
    A[4,1]=ks[1]*(35/(35+190+30))
    A[5,1]=ks[1]*(30/(35+190+30))  
    
    At=BoundLinDecompOp(
           map=function(t){ fX(t)*A },
           t_start,
           t_stop
           ) 
    
    Fc=BoundFc(inputFc,lag=lag,format="Delta14C")
    
    #mod=GeneralModel_14(t,
    mod=Model_14(t,
      At,
      ivList=C0,
      initialValF=ConstFc(F0_Delta14C,"Delta14C"),
      inputFluxes=inputFluxes,
      inputFc=Fc,
      c14DecayRate=lambda,
      #di=lambda,
      pass=pass
    )
    ### A Model Object that can be further queried 
    ##seealso<< \code{\link{ThreepParallelModel14}}, \code{\link{ThreepFeedbackModel14}} 
  }
  ,
  ex=function(){
    
    years=seq(1901,2010,by=0.5)
    
    Ex=GaudinskiModel14(
      t=years,
      ks=c(kr=1/3, koi=1/1.5, koeal=1/4, koeah=1/80, kA1=1/3, kA2=1/75, kM=1/110),
      inputFc=C14Atm_NH
    )
    R14m=getF14R(Ex)
    C14m=getF14C(Ex)
        
    plot(
      C14Atm_NH,
      type="l",
      xlab="Year",
      ylab=expression(paste(Delta^14,"C ","(\u2030)")),
      xlim=c(1940,2010)
    ) 
    lines(years,C14m,col=4)
    points(HarvardForest14CO2[1:11,1],HarvardForest14CO2[1:11,2],pch=19,cex=0.5)
    points(HarvardForest14CO2[12:173,1],HarvardForest14CO2[12:173,2],pch=19,col=2,cex=0.5)
    points(HarvardForest14CO2[158,1],HarvardForest14CO2[158,2],pch=19,cex=0.5)
    lines(years,R14m,col=2)
    legend(
      "topright",
      c("Delta 14C Atmosphere",
         "Delta 14C SOM", 
         "Delta 14C Respired"
      ),
      lty=c(1,1,1), 
      col=c(1,4,2),
      bty="n"
    )
    ## We now show how to bypass soilR s parameter sanity check if nacessary 
    ## (e.g in for parameter estimation ) in functions
    ## wchich might call it with unreasonable parameters
    years=seq(1800,2010,by=0.5)
    Ex=GaudinskiModel14(
      t=years,
      ks=c(kr=1/3,koi=1/1.5,koeal=1/4,koeah=1/80,kA1=1/3,kA2=1/75,kM=1/110),
      inputFc=C14Atm_NH,
      pass=TRUE
   )
  }
  )
