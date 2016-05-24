Longwave <-
function(emissivity,temp){
# daily longwave radiation based on the Sephan-Boltzman equation [kJ m-2 d-1]

#emissivity: [-]
#temp: temperature of the emitting body [C]

SBconstant<-0.00000490 #[kJ m-2 K-4 d-1]

tempK<-temp+273.15 #[degrees K]

return(emissivity*SBconstant*tempK^4)
}

