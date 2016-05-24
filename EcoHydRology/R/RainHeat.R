RainHeat <-
function(airtemp,rain){
# temperature added to the land from heat exchange with rain (usually in the context of snowmelt) [kJ m-2 d-1]

#airtemp: average dailiy air temperature [C]
#rain: depth of rainfall [m]

heatcapacity<-4.2  #heat capacity of water [kJ kg-1 C-1]
waterdensity<-1000 #density of water [kg m-3]

return(heatcapacity*waterdensity*rain*(airtemp-0))
}

