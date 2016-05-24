## =============================================================================
##
## The Daphnia model from Soetaert and Herman, 2009.
## a practical guide to ecological modelling,
## using R as a simulation platform. Springer.
## chapter 6
##
## implemented with 2 types of EVENTS:
##   transfer to new culture medium
##   moulting of the animals
##
## =============================================================================



library(deSolve)

#----------------------#
# the model equations: #
#----------------------#

model <- function(t, state, parameters)  {

  with(as.list(state), {  # unpack the state variables

    ## ingestion, size-dependent and food limited
    WeightFactor <- (IngestWeight - INDWEIGHT)/(IngestWeight - neonateWeight)
    MaxIngestion <- maxIngest * WeightFactor      # /day
    Ingestion    <- MaxIngestion * INDWEIGHT * FOOD / (FOOD + ksFood)

    Respiration  <- respirationRate * INDWEIGHT         # µgC/day
    Growth       <- Ingestion * assimilEff - Respiration

    ## Fraction of assimilate allocated to reproduction

    if (Growth <= 0 | INDWEIGHT < reproductiveWeight)
      Reproduction <- 0
    else {               # Fraction of growth allocated to reproduction.
      WeightRatio  <- reproductiveWeight/INDWEIGHT
      Reproduction <- maxReproduction * (1 - WeightRatio^2)
    }

    ## rate of change
    dINDWEIGHT <- (1 -Reproduction) * Growth
    dEGGWEIGHT <-      Reproduction * Growth
    dFOOD      <- -Ingestion * numberIndividuals

    ## the output, packed as a list
    list(c(dINDWEIGHT, dEGGWEIGHT, dFOOD),   # the rate of change
         c(Ingestion    = Ingestion,         # the ordinary output variables
           Respiration  = Respiration,
           Reproduction = Reproduction))
  })
}  # end of model

#---------------------------------------------------#
# Moulting weight loss and transfer to new culture  #
#---------------------------------------------------#

Eventfunc   <- function (t, state, parms) {
   with(as.list(state), {  # unpack the state variables
     if (t %in% MoultTime) {     # Moulting...
       ## Relationship moulting loss and length
       refLoss <- 0.24   #µgC
       cLoss   <- 3.1    #-

       ## Weight lost during molts depends allometrically on the organism length
       INDLength <- (INDWEIGHT /3.0)^(1/2.6)

       WeightLoss <- refLoss * INDLength^cLoss
       INDWEIGHT  <- INDWEIGHT - WeightLoss
       EGGWEIGHT  <- 0.
     }
     if (t %in% TransTime)   # New medium...
       FOOD <- foodInMedium

     return(c(INDWEIGHT, EGGWEIGHT, FOOD))
  })
}

#-----------------------#
# the model parameters: #
#-----------------------#

neonateWeight      <-  1.1    #µgC
reproductiveWeight <-  7.5    #µgC
maximumWeight      <- 60.0    #µgC

ksFood             <- 85.0    #µgC/l
IngestWeight       <-132.0    #µgC
maxIngest          <-  1.05   #/day
assimilEff         <-  0.8    #-

maxReproduction    <-  0.8    #-
respirationRate    <-  0.25   #/day

## Dilution parameters !
transferTime       <-    2    # Days
foodInMedium       <-  509    # µgC/l

instarDuration     <-  3.0    # days
numberIndividuals  <-   32    #   -

#-------------------------#
# the initial conditions: #
#-------------------------#

state <- c(
  INDWEIGHT = neonateWeight,       # µgC
  EGGWEIGHT = 0,                   # µgC    ! Total egg mass in a stage
  FOOD      = foodInMedium         # µgC
)

#----------------------#
# RUNNING the model:   #
#----------------------#
TimeEnd      <- 40                          # duration of simulation, days

times <- seq(0, TimeEnd, 0.1)               # output array

## when events are happening...
MoultTime <- seq(from = instarDuration, to = TimeEnd, by = instarDuration)
TransTime <- seq(from = transferTime, to = TimeEnd, by = transferTime)
EventTime <- sort(unique(c(MoultTime, TransTime)))

out <- ode(times = times, func = model, parms = NULL, y = state,
   events = list(func = Eventfunc, time = EventTime))

par(mfrow = c(2, 2),  oma = c(0, 0, 3, 0))   # set number of plots (mfrow) and margin size (oma)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot (out, which = c("FOOD", "INDWEIGHT", "EGGWEIGHT", "Ingestion"), type = "l",
      xlab = "time, days", ylab = c("gC/m3", "µgC", "µgC", "µgC/day"))
#main = "Food"              ,
#plot (out, which = , type = "l", main = "individual weight" , xlab = "time, days", ylab=)
#plot (out, which = , type = "l", main = "egg weight"        , xlab = "time, days", ylab=)
#plot (out, which = , type = "l", main = "Ingestion"         , xlab = "time, days", ylab=)

mtext(outer = TRUE, side = 3, "DAPHNIA model", cex = 1.5)

