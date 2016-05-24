par(ask=TRUE)

######################################################################################
# Parameters
######################################################################################
parameters <- list(             
                   t          = 15    , # degrees C
                   S          = 25    , # psu       
                   
                   k          = 0.4       , # 1/d	    proportionality factor for air-water exchange
                   rOx        = 0.0000003 , # mol-N/(kg*d)  maximal rate of oxic mineralisation
                   rNitri     = 0.0000002 , # mol-N/(kg*d)  maximal rate of nitrification 
                   rPP        = 0.0000006 , # mol-N/(kg*d)  maximal rate of primary production
                   
                   ksSumNH4   = 0.000001  , # mol-N/kg
                  
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
boxmodel <- function(timestep, currentstate, parameters)
{
  with (
        as.list(c(currentstate,parameters)),
        {        
          ae <- aquaenv(S=S, t=t, SumCO2=SumCO2, SumNH4=SumNH4, TA=TA, dsa=TRUE)
                                    
          ECO2    <- k * (ae$CO2_sat - ae$CO2)            
          EO2     <- k * (ae$O2_sat  - O2)                    
        
          RNit      <- rNitri 
          ROx       <- rOx 
        
          if ((timestep > a) && (timestep < b))
            {
              RPP <- rPP * (SumNH4/(ksSumNH4 + SumNH4))
            }
          else
            {
              RPP <- 0
            }
          
          dO2     <- EO2 - C_Nratio*ROx - 2*RNit + C_Nratio*RPP
          dNO3    <- RNit
          
          dSumCO2 <- ECO2 + C_Nratio*ROx - C_Nratio*RPP
          dSumNH4 <- ROx  - RNit - RPP

          dTA     <- ROx - 2*RNit - RPP

          
          # The DSA pH
          dH    <- (dTA - (dSumCO2*ae$dTAdSumCO2 + dSumNH4*ae$dTAdSumNH4))/ae$dTAdH
          DSApH <- -log10(H)

          # The DSA pH using pH dependent fractional stoichiometry (= using partitioning coefficients)
          rhoHECO2 <- ae$c2 + 2*ae$c3
          rhoHRNit <- 1 + ae$n1
          rhoHROx  <- C_Nratio * (ae$c2 + 2*ae$c3) - ae$n1
          rhoHRPP  <- -(C_Nratio * (ae$c2 + 2*ae$c3)) + ae$n1
          
          dH_ECO2  <- rhoHECO2*ECO2/(-ae$dTAdH)
          dH_RNit  <- rhoHRNit*RNit/(-ae$dTAdH)
          dH_ROx   <- rhoHROx*ROx  /(-ae$dTAdH)
          dH_RPP   <- rhoHRPP*RPP  /(-ae$dTAdH)

          dH_stoich   <- dH_ECO2 + dH_RNit + dH_ROx + dH_RPP
          DSAstoichpH <- -log10(H_stoich)       

          ratesofchanges <- c(dO2, dNO3, dSumNH4, dSumCO2, dTA, dH, dH_stoich)
          processrates   <- c(ECO2=ECO2, EO2=EO2, RNit=RNit, ROx=ROx, RPP=RPP)
          DSA            <- c(DSApH=DSApH, rhoHECO2=rhoHECO2, rhoHRNit=rhoHRNit, rhoHROx=rhoHROx,
                              rhoHRPP=rhoHRPP, dH_ECO2=dH_ECO2, dH_RNit=dH_RNit, dH_ROx=dH_ROx, dH_RPP=dH_RPP, DSAstoichpH=DSAstoichpH)
          
          return(list(ratesofchanges, list(processrates, DSA, ae)))
        }
        )
}



######################################################################################
# The model solution
######################################################################################
with (as.list(parameters),
      {
        H_init       <<- 10^(-(aquaenv(S=S, t=t, SumCO2=SumCO2_io, SumNH4=SumNH4_io, TA=TA_io, speciation=FALSE)$pH))
        initialstate <<- c(O2=O2_io, NO3=NO3_io, SumNH4=SumNH4_io, SumCO2=SumCO2_io, TA=TA_io, H=H_init, H_stoich=H_init)
        times        <<- c(0:modeltime)
        output       <<- as.data.frame(vode(initialstate, times, boxmodel, parameters, hmax=1))[-1,]        
      })


######################################################################################
# Output
######################################################################################

# plot the most important modelresults
what <- c("SumCO2", "TA", "SumNH4", "NO3", "ECO2", "EO2", "RNit", "ROx", "RPP", "dTAdH", "dTAdSumCO2", "dTAdSumNH4", "c1", "c2", "c3", "n1", "n2",
          "rhoHECO2", "rhoHRNit", "rhoHROx", "rhoHRPP", "dH_ECO2", "dH_RNit", "dH_ROx", "dH_RPP",
          "pH", "DSApH", "DSAstoichpH")
plot(aquaenv(ae=output, from.data.frame=TRUE), xval=output$time, what=what,  xlab="time/d", mfrow=c(6,5), size=c(20,13), newdevice=FALSE) 

# cumulative plot of the influences of the processes on the pH
par(mfrow=c(1,2))
what <- c("dH_ECO2", "dH_RNit", "dH_ROx", "dH_RPP")
plot(aquaenv(ae=output, from.data.frame=TRUE), xval=output$time, what=what, xlab="time/d", size=c(7,5), ylab="mol-H/(kg-soln*d)", legendposition="bottomright", cumulative=TRUE, newdevice=FALSE) 
 
# check if all three methods of calculating the pH yield the same result
ylim <- range(output$DSApH, output$DSAstoichpH, output$pH)
plot(output$DSApH, ylim=ylim, type="l")
par(new=TRUE)
plot(output$DSApH, ylim=ylim, type="l", col="red")
par(new=TRUE)
plot(output$DSAstoichpH, ylim=ylim, type="l", col="blue")




