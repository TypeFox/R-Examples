search.model.LM <- function(version = c("basic","latent","manifest"), kv, ..., nrep=2, tol1 = 10^-5, tol2 = 10^-10,out_se=FALSE){

# function that search for the global maximum of the log-likelihood
# vector of kv to try for
# nrep = number repetitions with random starting values
# version = model to be estimated ("basic" = basic LM model (est_lm_basic function); "latent" = LM model with covariates in the distribution of the latent process (est_lm_cov_latent function); "manifest" = LM model with covariates in the measurement model (est_lm_cov_maifest function))


out = vector("list",max(kv))
lkv = aicv = bicv = errv = rep(NA,max(kv))
#if(min(kv)>1) kv = c(1,kv)

version = match.arg(version)
for(k in kv){
	cat("***************************************************************************\n")
  	cat(k,"\n")
  	if(version=="basic") out[[k]] = try(est_lm_basic(...,k=k,start=0,tol=tol1))
	if(version=="latent") out[[k]] = try(est_lm_cov_latent(...,k=k,start=0,tol=tol1))	
	if(version=="manifest") out[[k]] = try(est_lm_cov_manifest(...,k=k,start=0,tol=tol1))
	
	if(!inherits(out[[k]],"try-error")){
  		errv[[k]] = FALSE
  	}else{
  		errv[[k]] = TRUE
  		if(k>1) out[[k]] = out[[k-1]]
 	}
 	lktrace = out[[k]]$lk
  	lkv[k] = out[[k]]$lk
  	aicv[k] = out[[k]]$aic
  	bicv[k] = out[[k]]$bic
	cat("lktrace = ",sort(lktrace),"\n")
 	cat("lk = ",lkv,"\n")
  	cat("aic = ",aicv,"\n")
  	cat("bic = ",bicv,"\n")
  	save(file="search.LM.temp.RData",out,aicv,bicv,errv,lkv)
  	if(k>1){ 
  		if(nrep==0){
  			cat("***************************************************************************\n")
  				cat(c(k,1),"\n")
    			if(version=="basic") outh = try(est_lm_basic(...,k=k,start=1,tol=tol1))	
    			if(version=="latent") outh = try(est_lm_cov_latent(...,k=k,start=1,tol=tol1))
    			if(version=="manifest") outh = try(est_lm_cov_manifest(...,k=k,start=1,tol=tol1))
    			if(!inherits(outh,"try-error")){
	    			lktrace = c(lktrace,outh$lk)
					if(outh$lk>out[[k]]$lk) out[[k]] = outh	
    			}  	
				lkv[k] = out[[k]]$lk
  				aicv[k] = out[[k]]$aic
    			bicv[k] = out[[k]]$bic
   		
   				cat("lktrace = ",sort(lktrace),"\n")
    			cat("lk = ",lkv,"\n")
    			cat("aic = ",aicv,"\n")
    			cat("bic = ",bicv,"\n")
    			save(file="search.LM.temp.RData",out,aicv,bicv,errv,lkv) 
  		}else{
  			for(h in 1:(nrep*(k-1))){
    			cat("***************************************************************************\n")
  				cat(c(k,h),"\n")
    			if(version=="basic") outh = try(est_lm_basic(...,k=k,start=1,tol=tol1))	
    			if(version=="latent") outh = try(est_lm_cov_latent(...,k=k,start=1,tol=tol1))
    			if(version=="manifest") outh = try(est_lm_cov_manifest(...,k=k,start=1,tol=tol1))
    			if(!inherits(outh,"try-error")){
	    			lktrace = c(lktrace,outh$lk)
					if(outh$lk>out[[k]]$lk) out[[k]] = outh	
    			}  	
				lkv[k] = out[[k]]$lk
  				aicv[k] = out[[k]]$aic
    			bicv[k] = out[[k]]$bic
   		
   				cat("lktrace = ",sort(lktrace),"\n")
    			cat("lk = ",lkv,"\n")
    			cat("aic = ",aicv,"\n")
    			cat("bic = ",bicv,"\n")
    			save(file="search.LM.temp.RData",out,aicv,bicv,errv,lkv) 		
  			}
  		}
  		if(version=="basic") outn = try(est_lm_basic(...,k=k,start=2,tol=tol2,piv=out[[k]]$piv,Pi=out[[k]]$Pi,Psi=out[[k]]$Psi,out_se=out_se))	
  	
  		if(version=="latent") outn = try(est_lm_cov_latent(...,k=k,start=2,tol=tol2,Psi=out[[k]]$Psi,Be=out[[k]]$Be,Ga=out[[k]]$Ga,out_se=out_se))   	
  	
  		if(version=="manifest") outn = try(est_lm_cov_manifest(...,k=k,start=2,tol=tol2,mu=out[[k]]$mu,al=out[[k]]$al,be=out[[k]]$be,la=out[[k]]$la,PI=out[[k]]$PI,rho=out[[k]]$rho,si=out[[k]]$si,out_se=out_se))
  		if(!inherits(outn,"try-error")){
  			lktrace = c(lktrace,outn$lk)
  			out[[k]] = outn		
  			out[[k]]$lktrace = lktrace
  			lkv[k] = out[[k]]$lk
  			aicv[k] = out[[k]]$aic
   	 		bicv[k] = out[[k]]$bic
  		}
  }	
  save(file="search.LM.temp.RData",out,aicv,bicv,errv,lkv)
}
	
out = list(out.single=out,aicv=aicv,bicv=bicv,lkv=lkv,errv=errv,kv=kv,call=match.call())
class(out)="LMsearch"	
return(out)
}