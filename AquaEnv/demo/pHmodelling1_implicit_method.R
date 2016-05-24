par(ask=TRUE)


what   <- c("SumCO2", "TA", "Rc", "Rp", "omega_calcite", "pH")
mfrow  <- c(2,3)
size   <- c(15,7)
           

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
          ae    <- aquaenv(S=S, t=t, SumCO2=SumCO2, TA=TA, SumSiOH4=0, SumBOH3=0, SumH2SO4=0, SumHF=0)      
          
          Rc    <- kc * ((ae$CO2_sat) - (ae$CO2)) 
          Rp    <- kp * (1-ae$omega_calcite)^n               

          dSumCO2 <- Rc - Rp
          dTA     <- -2*Rp
          
          ratesofchanges <- c(dSumCO2, dTA)
          
          processrates   <- c(Rc=Rc, Rp=Rp)
          
          return(list(ratesofchanges, list(processrates, ae)))
        }
        )
}


with (as.list(parameters),
      {
        initialstate <<- c(SumCO2=SumCO2_t0, TA=TA_t0)
        times        <<- seq(0,modeltime,(modeltime/outputsteps))       
        output       <<- as.data.frame(vode(initialstate,times,boxmodel,parameters, hmax=1))
      })



plot(aquaenv(ae=output, from.data.frame=TRUE), xval=output$time, xlab="", mfrow=mfrow, size=size, what=what, newdevice=FALSE) 




