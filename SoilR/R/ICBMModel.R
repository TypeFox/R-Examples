#
# vim:set ff=unix expandtab ts=2 sw=2:
ICBMModel<-structure(
    function #Implementation of the Introductory Carbon Balance Model (ICBM)
    ### This function is an implementation of the Introductory Carbon Balance Model (ICBM).
    ### This is simply a two pool model connected in series.
    ##references<< Andren, O. and T. Katterer. 1997. ICBM: The Introductory Carbon Balance Model 
    ##for Exploration of Soil Carbon Balances. Ecological Applications 7:1226-1236.

    (t, ##<< A vector containing the points in time where the solution is sougth.
     ks=c(k1=0.8,k2=0.00605), ##<< A vector of length 2 with the decomposition rates for the young and the old pool.
     h=0.13, ##<< Humufication coefficient (transfer rate from young to old pool).
     r=1.32, ##<< External (environmental or edaphic) factor. 
     c0=c(Y0=0.3,O0=3.96), ##<< A vector of length 2 with the initial value of carbon stocks in the young and old pool. 
     In=0,   ##<< Mean annual carbon input to the soil. 
     solver=deSolve.lsoda.wrapper, ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
     pass=FALSE  ##<< if TRUE forces the constructor to create the model even if it is invalid 
     )
    { 
      t_start=min(t)
      t_end=max(t)
     if(length(ks)!=2) stop("The vector of decomposition rates is not of length = 2")
     if(length(c0)!=2) stop("The vector with initial conditions is not of length = 2")
     
     A=diag(-ks)
     A[2,1]=ks[1]*h
     Ar=A*r
     inputFluxes=BoundInFlux(
        function(t){matrix(nrow=nrow(A),ncol=1,c(In,0))},
        t_start,
        t_end
     )
     Af=BoundLinDecompOp(map=function(t0){Ar},t_start,t_end)
     Mod=GeneralModel(t=t,A=Af,c0,inputFluxes,solver,pass)
     return(Mod)
 
     ##seealso<< \code{\link{TwopSeriesModel}} 
    }
    ,
    ex=function(){
        # This example reproduces the simulations 
        # presented in Table 1 of Andren and Katterer (1997).
        # First, the model is run for different values of the 
        # parameters representing different field experiments. 
        times=seq(0,20,by=0.1)
        Bare=ICBMModel(t=times) #Bare fallow
        pNpS=ICBMModel(t=times, h=0.125, r=1,    c0=c(0.3,4.11),  In=0.19+0.095) #+N +Straw
        mNpS=ICBMModel(t=times, h=0.125, r=1.22, c0=c(0.3, 4.05), In=0.19+0.058) #-N +Straw
        mNmS=ICBMModel(t=times, h=0.125, r=1.17, c0=c(0.3, 3.99), In=0.057) #-N -Straw
        pNmS=ICBMModel(t=times, h=0.125, r=1.07, c0=c(0.3, 4.02), In=0.091) #+N -Straw
        FM=ICBMModel(t=times, h=0.250, r=1.10, c0=c(0.3, 3.99), In=0.19+0.082) #Manure
        SwS=ICBMModel(t=times, h=0.340, r=0.97, c0=c(0.3, 4.14), In=0.19+0.106) #Sewage Sludge
        SS=ICBMModel(t=times, h=0.125, r=1.00, c0=c(0.25, 4.16), In=0.2)  #Steady State

        #The amount of carbon for each simulation is recovered with the function getC
        CtBare=getC(Bare)
        CtpNpS=getC(pNpS)
        CtmNpS=getC(mNpS)
        CtmNmS=getC(mNmS)
        CtpNmS=getC(pNmS)
        CtFM=getC(FM)
        CtSwS=getC(SwS)
        CtSS=getC(SS)

        #This plot reproduces Figure 1 in Andren and Katterer (1997)
        plot(times,
          rowSums(CtBare),
          type="l",
          ylim=c(0,8),
          xlim=c(0,20),
          ylab="Topsoil carbon mass (kg m-2)",
          xlab="Time (years)"
        )
        lines(times,rowSums(CtpNpS),lty=2)
        lines(times,rowSums(CtmNpS),lty=3)
        lines(times,rowSums(CtmNmS),lty=4)
        lines(times,rowSums(CtpNmS),lwd=2)
        lines(times,rowSums(CtFM),lty=2,lwd=2)
        lines(times,rowSums(CtSwS),lty=3,lwd=2)
        #lines(times,rowSums(CtSS),lty=4,lwd=2)
        legend("topleft",
          c("Bare fallow",
            "+N +Straw",
            "-N +Straw",
            "-N -Straw",
            "+N -Straw",
            "Manure",
           "Sludge"
          ),
          lty=c(1,2,3,4,1,2,3),
          lwd=c(1,1,1,1,2,2,2),
          bty="n"
        )
 
}
)
