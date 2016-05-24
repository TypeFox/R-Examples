working.compvariogmodels1 <-
function(gammadat,models,range0)
{
 ## gammadat : x : distance, gamma empirical semivariance
 ## models: variogram models to fit
 ## range0 : first peak or range
 ## value : better fitting models weighted least square
 ## note to Dong July 29: add nugget effects
 lsout <- vector('list',length(models))
 lsfit <- rep(NA,length(models))
 coefs <- matrix(NA,length(models),2)
 for(i in 1:length(models)){
	##cat(models[i],'\n')
	f <- get(paste('v',models[i],sep=''))
	xx <- with(gammadat,f(x,1,range0,0))
	## catch <- try(lm(gamma ~ xx-1,data=gammadat,weight=n),silent=T)
	catch <- try(lm(gamma ~ xx,data=gammadat,weights=gammadat$n),silent=T)
    chk <- inherits(catch,'try-error')
	if(!chk) {
		lsout[[i]]<- catch
		lsfit[i] <- deviance(lsout[[i]])
		coefs[i,] <- coef(catch)
	}
 }
 tmp <- apply(coefs,1,function(v) 
	{ out <- all(!is.na(v))
	  if(out) {
		out <- all(v>0)
	  }
	  out }
     )
 vii <- which(tmp)
 ii <- vii[which.min(lsfit[tmp])]
 if(length(ii)==0) {
	ii <- 1  ## chose exp by default
	xx <- vexpn(gammadat$x,1,range0,0)
	lsout[[ii]] <- lm(gamma ~ xx -1, data=gammadat,weights=gammadat$n)
	coefs[ii,] <- c(0,coef(lsout[[ii]]))
 }
 list(parms=coefs[ii,],model=models[ii])
}
