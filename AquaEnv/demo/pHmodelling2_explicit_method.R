par(ask=TRUE)


what   <- c("SumCO2", "TA", "Rc", "Rp", "omega_calcite", "pH", "dHRc", "dHRp")
mfrow  <- c(3,3)
size   <- c(15,10)
           

parameters <- list(             
                   S           = 35        , # psu       
                   t           = 15        , # degrees C

                   SumCO2_t0   = 0.002     , # mol/kg-soln  (comparable to Wang2005)
                   TA_t0       = 0.0022    , # mol/kg-soln  (comparable to Millero1998)

                   kc          = 0.5       , # 1/d	         proportionality factor for air-water exchange
                   kp          = 0.000001  , # mol/(kg-soln*d)	 max rate of calcium carbonate precipitation
                   n           = 2.0       , # -                 exponent for kinetic rate law of precipitation
                                      
                   modeltime   = 20        , # d              duration of the model
                   outputsteps = 100         #                number of outputsteps
                   )


boxmodel <- function(timestep, currentstate, parameters)
{
  with (
        as.list(c(currentstate,parameters)),
        {        
          ae    <- aquaenv(S=S, t=t, SumCO2=SumCO2, pH=-log10(H), SumSiOH4=0, SumBOH3=0, SumH2SO4=0, SumHF=0, dsa=TRUE)
                   
          Rc    <- kc * ((ae$CO2_sat) - (ae$CO2)) 
          Rp    <- kp * (1-ae$omega_calcite)^n               

          dSumCO2 <- Rc - Rp

          dHRc    <- (      -(ae$dTAdSumCO2*Rc   ))/ae$dTAdH
          dHRp    <- (-2*Rp -(ae$dTAdSumCO2*(-Rp)))/ae$dTAdH
          dH      <- dHRc + dHRp
          
          ratesofchanges <- c(dSumCO2, dH)
          
          processrates   <- c(Rc=Rc, Rp=Rp)
          outputvars     <- c(dHRc=dHRc, dHRp=dHRp)
          
          return(list(ratesofchanges, list(processrates, outputvars, ae)))
        }
        )
}


with (as.list(parameters),
      {
        aetmp <- aquaenv(t=t, S=S, SumCO2=SumCO2_t0, TA=TA_t0, SumSiOH4=0, SumBOH3=0, SumH2SO4=0, SumHF=0)
        H_t0  <- 10^(-aetmp$pH)
        
        initialstate <<- c(SumCO2=SumCO2_t0, H=H_t0)
        times        <<- seq(0,modeltime,(modeltime/outputsteps))       
        output       <<- as.data.frame(vode(initialstate,times,boxmodel,parameters, hmax=1))
      })


plot(aquaenv(ae=output, from.data.frame=TRUE), xval=output$time, xlab="time/d", mfrow=mfrow, size=size, what=what, newdevice=FALSE) 



