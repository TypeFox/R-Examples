SatVaporPressure <-
function(T_C){  
# saturated vapor pressure at a given temperature (kPa)
#T_C: temperature [C]
return(0.611 * exp((17.3*T_C)/(237.2+T_C)))
}
