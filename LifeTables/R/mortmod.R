mortmod <-
function(pattern, sex="female", alpha=0){
		x <- seq(0,5,.5)
	 	csd.prop <- 1-pexp(x, rate=.75) 
 		f.csd.weight <- approxfun(x, csd.prop)
 		
 		w.ave <- function(csd, oad, csd.weight){
		oad.weight <- 1-csd.weight
		dev.out <- (csd*csd.weight) + (oad*oad.weight)
		return(dev.out)
			}

		if(alpha < 0){
			to.subtract <- w.ave(csd=get("lo.devs",envir=.GlobalEnv)[pattern,], oad=get("lo.devs",envir=.GlobalEnv)[nrow(get("lo.devs",envir=.GlobalEnv)),], 
				csd.weight=f.csd.weight(abs(alpha)))
			model.patt <- get("averages.smooth",envir=.GlobalEnv)[,pattern] + alpha*to.subtract} else {
			to.add <- w.ave(csd=get("hi.devs",envir=.GlobalEnv)[pattern,], oad=get("hi.devs",envir=.GlobalEnv)[nrow(get("hi.devs",envir=.GlobalEnv)),], 
				csd.weight=f.csd.weight(abs(alpha)))
			model.patt <- get("averages.smooth",envir=.GlobalEnv)[,pattern] + alpha*to.add
			}
			
			model.patt[model.patt>0] <- -.0001
			if(sex=="male"){
			return(structure(model.patt[1:24], class='MLT'))
				}
				
			if(sex=="female"){
			return(structure(model.patt[25:48], class='MLT'))
				}
		}

