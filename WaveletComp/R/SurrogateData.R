SurrogateData <-
function(x, method = "white.noise", 
                          params=list(AR     = list(p=1),
                                      ARIMA  = list(p=1, q=1, include.mean=T, sd.fac=1, trim = F, trim.prop = 0.01)
#                                       ,
#                                       meboot = list(trim = 0.1, force.clt = F, expand.sd = T, fiv = 5)
                                     ) 
                          ){
                          
  if(method == "white.noise")  x.sur <- rnorm(length(x)) 
  if(method == "shuffle")      x.sur <- sample(x, length(x)) 
  if(method == "Fourier.rand") x.sur <- FourierRand(x) 
  
  if(method == "AR")           { 
 
     x.sur <- AR(x, params = params) 
     
  } 
  
#   if(method == "meboot")       { 
#   
#      trim      = params$meboot$trim
#      force.clt = params$meboot$force.clt
#      expand.sd = params$meboot$expand.sd
#      fiv       = params$meboot$fiv
#      
#      x.sur <- meboot(x, reps=2, trim = trim, force.clt = force.clt, expand.sd = expand.sd, fiv = fiv)$ensemble[,1]
#      
#   }
  
  if(method == "ARIMA")         {
  
     x.sur <- ARIMA(x, params = params)
 
  }
  
  return(invisible(x.sur))
}
