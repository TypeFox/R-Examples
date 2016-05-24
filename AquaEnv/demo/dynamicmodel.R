par(ask=TRUE)

######################################################################################
# Parameters
######################################################################################
parameters <- list(             
                   S          = 25    , # psu       
                   t_min      = 5     , # degrees C
                   t_max      = 25    , # degrees C
                   d          = 10    , # m
                   
                   k          = 0.4       , # 1/d	      proportionality factor for air-water exchange
                   rOx        = 0.0000003 , # mol-N/(kg*d)  maximal rate of oxic mineralisation
                   rNitri     = 0.0000002 , # mol-N/(kg*d)  maximal rate of nitrification 
                   rPP        = 0.000006  , # mol-N/(kg*d)  maximal rate of primary production
                   
                   ksDINPP    = 0.000001  , # mol-N/kg
                   ksNH4PP    = 0.000001  , # mol-N/kg
                   
                   D          = 0.1       , #   1/d            (dispersive) transport coefficient
                   
                   O2_io      = 0.000296  , # mol/kg-soln 
                   NO3_io     = 0.000035  , # mol/kg-soln 
                   SumNH4_io  = 0.000008  , # mol/kg-soln 
                   SumCO2_io  = 0.002320  , # mol/kg-soln 
                   TA_io      = 0.002435  , # mol/kg-soln 
                   
                   C_Nratio     = 8       , # mol C/mol N     C:N ratio of organic matter
                   
                   a           = 30       , # timestep from which PP begins     
                   b           = 50       , # timestep where PP shuts off again
                   
                   modeltime   = 100        # duration of the model
                   )


######################################################################################
# The model function
######################################################################################
Waddenzeebox <- function(timestep, currentstate, parameters)
{
  with (
        as.list(c(currentstate,parameters)),
        {
          t <- c(seq(t_min, t_max, (t_max-t_min)/(modeltime/2)), seq(t_max, t_min, -(t_max-t_min)/(modeltime/2)))[[round(timestep)+1]]
          
          ae <- aquaenv(S=S, t=t, SumCO2=SumCO2, SumNH4=SumNH4, TA=TA)
                                    
          ECO2    <- k * (ae$CO2_sat - ae$CO2)            
          EO2     <- k * (ae$O2_sat  - O2)             
         
          TO2     <- D*(O2_io     - O2)
          TNO3    <- D*(NO3_io    - NO3)
          TSumNH4 <- D*(SumNH4_io - SumNH4)
          TTA     <- D*(TA_io     - TA)
          TSumCO2 <- D*(SumCO2_io - SumCO2)
          
          RNit      <- rNitri 

          ROx       <- rOx 
          ROxCarbon <- ROx * C_Nratio

          pNH4PP <- 0
          RPP <- 0
          
          if ((timestep > a) && (timestep < b))
            {
              RPP    <- rPP * ((SumNH4+NO3)/(ksDINPP + (SumNH4+NO3)))
              pNH4PP <- 1 - (ksNH4PP/(ksNH4PP + SumNH4))
            }
          else
            {
              RPP <- 0
            }
          RPPCarbon <- RPP * C_Nratio
          
          dO2     <- TO2     + EO2 - ROxCarbon - 2*RNit  + (2-2*pNH4PP)*RPP + RPPCarbon
          dNO3    <- TNO3    + RNit -(1-pNH4PP)*RPP

          dSumCO2 <- TSumCO2 + ECO2 + ROxCarbon - RPPCarbon
          dSumNH4 <- TSumNH4 + ROx  - RNit - pNH4PP*RPP
                    
          dTA     <- TTA     + ROx - 2*RNit -(2*pNH4PP-1)*RPP 

          ratesofchanges <- c(dO2, dNO3, dSumNH4, dSumCO2, dTA)
          transport      <- c(TO2=TO2, TNO3=TNO3, TSumNH4=TSumNH4, TTA=TTA, TSumCO2=TSumCO2)
          airseaexchange <- c(ECO2=ECO2, EO2=EO2)
          
          return(list(ratesofchanges, list(transport, airseaexchange, ae)))
        }
        )
}



######################################################################################
# The model solution
######################################################################################
with (as.list(parameters),
      {
        initialstate <<- c(O2=O2_io, NO3=NO3_io, SumNH4=SumNH4_io, SumCO2=SumCO2_io, TA=TA_io)
        times        <<- c(0:modeltime)
        output       <<- as.data.frame(vode(initialstate,times,Waddenzeebox,parameters, hmax=1))[-1,]        
      })


######################################################################################
# Output
######################################################################################
plot(aquaenv(ae=output, from.data.frame=TRUE), xval=output$time, xlab="time/d", mfrow=c(10,11), newdevice=FALSE) 



