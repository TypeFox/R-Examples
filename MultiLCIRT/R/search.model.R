search.model <- function(S, yv = rep(1,ns), kv, X = NULL, link = 0, disc = 0,
	difl = 0, multi = 1:J, fort = FALSE, tol = 10^-10, nrep = 2, glob = FALSE, disp=FALSE){

# function that search for the global maximum of the log-likelihood
# vector of kv to try for
# nrep = number repetitions with random starting values
ns = dim(S)[1]
J = dim(S)[2]
out = vector("list",max(kv))
lkv = aicv = bicv = entv = necv = errv = rep(NA,max(kv))
if(min(kv)>1) kv = c(1,kv)
for(k in kv){
  cat("***************************************************************************\n")
  cat(k,"\n")
  out[[k]] = try(est_multi_poly(S=S,yv=yv,k=k,X=X,start=0,link=link,disc=disc,difl=difl,multi=multi,fort=fort,tol=tol,glob=glob,disp=disp))
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
  entv[k] = out[[k]]$ent
  if(k==1) necv[1] = 1
  else if(1%in%kv) necv[k] = entv[k]/(lkv[k]-lkv[1])
  cat("lktrace = ",sort(lktrace),"\n")
  cat("lk = ",lkv,"\n")
  cat("aic = ",aicv,"\n")
  cat("bic = ",bicv,"\n")
  cat("ent = ",entv,"\n")
  cat("nec = ",necv,"\n")
  save(file="search.model.temp.RData",out,aicv,bicv,entv,necv,errv,lkv)
  if(k>1) for(h in 1:(nrep*(k-1))){
    cat("***************************************************************************\n")
  	cat(c(k,h),"\n")
    outh = try(est_multi_poly(S=S,yv=yv,k=k,X=X,start=1,link=link,disc=disc,difl=difl,multi=multi,fort=fort,tol=tol,glob=glob,disp=disp))
    if(!inherits(outh,"try-error")){
	    lktrace = c(lktrace,outh$lk)
		if(outh$lk>out[[k]]$lk) out[[k]] = outh	
    }  	
	lkv[k] = out[[k]]$lk
  	aicv[k] = out[[k]]$aic
    bicv[k] = out[[k]]$bic
    entv[k] = out[[k]]$ent
    if(1%in%kv) necv[k] = entv[k]/(lkv[k]-lkv[1])
    cat("lktrace = ",sort(lktrace),"\n")
    cat("lk = ",lkv,"\n")
    cat("aic = ",aicv,"\n")
    cat("bic = ",bicv,"\n")
    cat("ent = ",entv,"\n")
    cat("nec = ",necv,"\n")
	save(file="search.model.temp.RData",out,aicv,bicv,entv,necv,errv,lkv)
  }
  out[[k]]$lktrace = lktrace
  save(file="search.model.temp.RData",out,aicv,bicv,entv,necv,errv,lkv)
}
out = list(out.single=out,aicv=aicv,bicv=bicv,entv=entv,necv=necv,lkv=lkv,errv=errv)

}