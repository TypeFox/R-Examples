ltsbaseSummary <-
function(object){

res=list(mses=object$list.mse, par=object$list.bias.par, coef=object$list.coef.all)
if(which.min(res$mses)==6) {cat("biasing parameter at best mse is"); print(res$par[4]);cat("corresponding coefficients", sep="\n");print(res$coef[,6])} else 
	if(which.min(res$mses)==5) {cat("biasing parameter at best mse is"); print(res$par[3]);cat("corresponding coefficients", sep="\n");print(res$coef[,5])} else 
		if(which.min(res$mses)==3) {cat("biasing parameter at best mse is"); print(res$par[2]);cat("corresponding coefficients", sep="\n");print(res$coef[,3])} else 
			if(which.min(res$mses)==2) {cat("biasing parameter at best mse is"); print(res$par[1]);cat("corresponding coefficients", sep="\n");print(res$coef[,2])} 
				if(which.min(res$mses)==1 || which.min(res$mses)==4) {print("No biasing parameters for OLS or LTS") } 

cat("best mse ", sep="\n")
best.mse=min(res$mses); print(best.mse)	
}
