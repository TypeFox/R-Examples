
print_uncertainty <- function(model,T,type="pn",...){
	
	d <- model@d
	if(d==1){
    return(print_uncertainty_1d(model=model,T=T,type=type,...))
	}else if(d==2){
    return(print_uncertainty_2d(model=model,T=T,type=type,...))
	}else{
    return(print_uncertainty_nd(model=model,T=T,type=type,...))
	}
  
}
