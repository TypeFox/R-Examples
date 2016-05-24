`growth2vgc.fnc` <-
function(growth) {
  require(zipfR, quietly = TRUE)
  return(vgc(N  = growth@data$data$Tokens, 
             V  = growth@data$data$Types, 
             Vm = list(growth@data$data$HapaxLegomena,
                        growth@data$data$DisLegomena, 
                        growth@data$data$TrisLegomena)))

}

