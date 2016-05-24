

searpas <- function(vw,step,b,delta,funcpa,res.out.error){

	goto50  <- function(step,vlw2,fi1,fi2,fi3,b,delta,funcpa){
		vm <- vlw2-(step*(fi1-fi3))/(2*(fi1-2*fi2+fi3)) 
		fim <- valfpa(vm,b,delta,funcpa)
		return(list(vm=vm,fim=fim))
	}
	vlw1 <- log(vw)
	vlw2 <- vlw1+step
	fi1 <- valfpa(vlw1,b,delta,funcpa)
	fi2 <- valfpa(vlw2,b,delta,funcpa)
	
	if((sum(!is.finite(fi1)) > 0) || (sum(!is.finite(fi2)) > 0)){
		cat("Probably too much accuracy requested...\n")
		cat("Last step values :\n")
		cat("      b :",res.out.error$old.b,"\n")
		cat("      likelihood :",res.out.error$old.rl,"\n")
		cat("      Convergence criteria: parameters stability=", res.out.error$old.ca, "\n")
		cat("                          : likelihood stability=", res.out.error$old.cb, "\n") 
		cat("                          : best relative distance to maximum obtained (RDM)=", res.out.error$old.dd, "\n")
		stop("")	
	}
	
	if((fi2 >= fi1)){
		vlw3 <- vlw2
		vlw2 <- vlw1
		fi3 <- fi2
		fi2 <- fi1
		step <- -step
		vlw1 <- vlw2+step
		fi1 <- valfpa(vlw1,b,delta,funcpa)
		gt50 <- goto50(step,vlw2,fi1,fi2,fi3,b,delta,funcpa)
		vm <- gt50$vm
		fim <- gt50$fim
		if(is.na(fim)) fim <- 10E10
		if(fim <= fi2){
			vw <- exp(vm)
		}else{
			vm <- vlw2
			fim <- fi2
			vw <- exp(vm)
		}

	}else{
		vlw <- vlw1
		vlw1 <- vlw2
		vlw2 <- vlw
		fim <- fi1
		fi1 <- fi2
		fi2 <- fim
		
		for(i in 1:40){
			vlw3 <- vlw2
			vlw2 <- vlw1
			fi3 <- fi2
			fi2 < fi1
			vlw1=vlw2+step
			fi1 <- valfpa(vlw1,b,delta,funcpa)
			if(fi1 > fi2){
				gt50 <- goto50(step,vlw2,fi1,fi2,fi3,b,delta,funcpa) 
				out <- 1
				break
			}
			if(fi1 == fi2){
				fim <- fi2
				vm <- vlw2
				vw <- exp(vm)
				out <- 1
				break
			}
		}
	}	
	return(list(vw=vw,fim=fim))
}