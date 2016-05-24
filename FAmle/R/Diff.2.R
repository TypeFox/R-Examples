Diff.2 <-
function(k,i,model,p,ln=FALSE)
{
	param <- model$par.hat
	param[i] <- k
	if(!ln) return(distr(p,model$dist,param,'q'))
	else return(log(distr(p,model$dist,param,'q')))
}

