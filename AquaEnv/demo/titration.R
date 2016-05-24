par(ask=TRUE)

#######################################################################################################
# Titration with HCl
#######################################################################################################
S  <- 35
t  <- 15

SumCO2 <- 0.003500
SumNH4 <- 0.000020

mass_sample  <- 0.01 # the mass of the sample solution in kg
mass_titrant <- 0.02 # the total mass of the added titrant solution in kg
conc_titrant <- 0.01 # the concentration of the titrant solution in mol/kg-soln
S_titrant    <- 0.5  # the salinity of the titrant solution (the salinity of a solution with a ionic strength of 0.01 according to: I = (19.924 S) / (1000 - 1.005 S)
steps        <- 50   # the amount of steps the mass of titrant is added in 
type         <- "HCl"

pHstart <- 11.3


ae <- titration(aquaenv(S=S, t=t, SumCO2=SumCO2, SumNH4=SumNH4, pH=pHstart), mass_sample, mass_titrant, conc_titrant, S_titrant, steps, type)


# plotting everything
plot(ae, xval=ae$delta_mass_titrant, xlab="HCl solution added [kg]", mfrow=c(10,10), newdevice=FALSE)


# plotting selectively
size  <- c(12,8) #inches
mfrow <- c(4,4)
what  <- c("TA", "pH", "CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F", "fCO2")

plot(ae, xval=ae$delta_mass_titrant, xlab="HCl solution added [kg]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)
plot(ae, xval=ae$pH, xlab="free scale pH", what=what, size=size, mfrow=mfrow, newdevice=FALSE)


# different x values
plot(ae, xval=ae$delta_conc_titrant, xlab="[HCl] offset added [mol/kg-soln]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)
plot(ae, xval=ae$delta_moles_titrant, xlab="HCl added [mol]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)



# bjerrum plots
par(mfrow=c(1,1))
plot(ae, bjerrum=TRUE, newdevice=FALSE)

what  <- c("CO2", "HCO3", "CO3")
plot(ae, what=what, bjerrum=TRUE, newdevice=FALSE)
plot(ae, what=what, bjerrum=TRUE, lwd=4, palette=c("cyan", "magenta", "yellow"), bg="gray", legendinset=0.1, legendposition="topleft", newdevice=FALSE)


what  <- c("CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F")
plot(ae, what=what, bjerrum=TRUE, log=TRUE, newdevice=FALSE)
plot(ae, what=what, bjerrum=TRUE, log=TRUE, ylim=c(-6,-1), legendinset=0, lwd=3, palette=c(1,3,4,5,6,colors()[seq(100,250,6)]), newdevice=FALSE)





#######################################################################################################
# Titration with NaOH
#######################################################################################################
S  <- 35
t  <- 15

SumCO2 <- 0.003500
SumNH4 <- 0.000020

mass_sample  <- 0.01 # the mass of the sample solution in kg
mass_titrant <- 0.02 # the total mass of the added titrant solution in kg
conc_titrant <- 0.01 # the concentration of the titrant solution in mol/kg-soln
S_titrant    <- 0.5  # the salinity of the titrant solution
steps        <- 50   # the amount of steps the mass of titrant is added in 
type         <- "NaOH"

pHstart <- 2

ae <- titration(aquaenv(S=S, t=t, SumCO2=SumCO2, SumNH4=SumNH4, pH=pHstart), mass_sample, mass_titrant, conc_titrant, S_titrant, steps, type)



# plottinge everything
plot(ae, xval=ae$delta_mass_titrant, xlab="NaOH solution added [kg]", mfrow=c(10,10), newdevice=FALSE)


# plotting selectively
size  <- c(12,8) #inches
mfrow <- c(4,4)
what  <- c("TA", "pH", "CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F", "fCO2")

plot(ae, xval=ae$delta_mass_titrant, xlab="NaOH solution added [kg]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)
plot(ae, xval=ae$pH, xlab="free scale pH", what=what, size=size, mfrow=mfrow, newdevice=FALSE)



# bjerrum plots
par(mfrow=c(1,1))
what  <- c("CO2", "HCO3", "CO3")
plot(ae, what=what, bjerrum=TRUE, newdevice=FALSE)

what  <- c("CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F")
plot(ae, what=what, bjerrum=TRUE, log=TRUE, ylim=c(-6,-1), legendinset=0, lwd=3, palette=c(1,3,4,5,6,colors()[seq(100,250,6)]), newdevice=FALSE)




#########################################################################################################################################################################
# titration with a titrant with very high concentrations and a very large sample volume: the volume and salinity corrections do not matter: looks like "simpletitration"
#########################################################################################################################################################################
S <- 35
t <- 15

SumCO2 <- 0.003500
SumNH4 <- 0.000020

mass_sample  <- 100  # the mass of the sample solution in kg
mass_titrant <- 0.5  # the total mass of the added titrant solution in kg
conc_titrant <- 3    # the concentration of the titrant solution in mol/kg-soln
S_titrant    <- 0.5  # the salinity of the titrant solution (the salinity of a solution with a ionic strength of 0.01 according to: I = (19.924 S) / (1000 - 1.005 S)
steps        <- 100  # the amount of steps the mass of titrant is added in 
type         <- "HCl"

pHstart <- 11.3


ae <- titration(aquaenv(S=S, t=t, SumCO2=SumCO2, SumNH4=SumNH4, pH=pHstart), mass_sample, mass_titrant, conc_titrant, S_titrant, steps, type)



# plotting everything
plot(ae, xval=ae$delta_mass_titrant, xlab="HCl solution added [kg]", mfrow=c(10,10), newdevice=FALSE)


# plotting selectively
size  <- c(12,8) #inches
mfrow <- c(4,4)
what  <- c("TA", "pH", "CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F", "fCO2")

plot(ae, xval=ae$delta_mass_titrant, xlab="HCl solution added [kg]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)
plot(ae, xval=ae$pH, xlab="free scale pH", what=what, size=size, mfrow=mfrow, newdevice=FALSE)


# different x values
plot(ae, xval=ae$delta_conc_titrant, xlab="[HCl] offset added [mol/kg-soln]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)
plot(ae, xval=ae$delta_moles_titrant, xlab="HCl added [mol]", what=what, size=size, mfrow=mfrow, newdevice=FALSE)



# bjerrum plots
par(mfrow=c(1,1))
plot(ae, bjerrum=TRUE, newdevice=FALSE)

what  <- c("CO2", "HCO3", "CO3")
plot(ae, what=what, bjerrum=TRUE, newdevice=FALSE)
plot(ae, what=what, bjerrum=TRUE, lwd=4, palette=c("cyan", "magenta", "yellow"), bg="gray", legendinset=0.1, legendposition="topleft", newdevice=FALSE)


what  <- c("CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F")
plot(ae, what=what, bjerrum=TRUE, log=TRUE, newdevice=FALSE)
plot(ae, what=what, bjerrum=TRUE, log=TRUE, ylim=c(-6,-1), legendinset=0, lwd=3, palette=c(1,3,4,5,6,colors()[seq(100,250,6)]), newdevice=FALSE)







