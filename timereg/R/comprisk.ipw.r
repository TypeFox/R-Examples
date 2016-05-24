comprisk.ipw <- function(compriskformula,glmformula,data=sys.parent(),cause=1,
			 max.clust=NULL,ipw.se=FALSE,...)
{
 ## {{{ 
  ggl <- glm(glmformula,family='binomial',data=data)
  mat <-  model.matrix(glmformula,data=data);
  glmcovs <- attr(ggl$terms,"term.labels")
  data$ppp <- predict(ggl,type='response')

  dcc <- data[ggl$y==1,]
  ppp <- dcc$ppp
  udca <- comp.risk(compriskformula,data=dcc,cause=cause,weights=1/ppp,n.sim=0,
		    max.clust=max.clust,...)  
  ### iid of beta for comprisk model 
  compriskiid <- udca$gamma.iid

if (ipw.se==TRUE)  { ## {{{ 
###requireNamespace("lava"); 
###requireNamespace("NumDeriv"); 
	glmiid <-   lava::iid(ggl)
	mat <- mat[ggl$y==1,]
	par <- coef(ggl)

	compriskalpha <- function(par)
	{ ## {{{ 
	  rr <- mat %*% par
	  pw <- c(exp(rr)/(1+exp(rr)))
	  assign("pw",pw,envir=environment(compriskformula))
	  ud <- comp.risk(compriskformula,data=dcc,cause=cause,
			  weights=1/pw,
			  est=udca$cum,gamma=udca$gamma,Nit=1,n.sim=0,...)  
	  ud$scores$gamscore
	} ## }}} 

	DU <-  numDeriv::jacobian(compriskalpha,par)
	IDU <-  udca$Dscore.gamma %*% DU 
	alphaiid <-t( IDU %*% t(glmiid))
	###
	iidfull <- alphaiid
	###
	iidfull[ggl$y==1,] <- compriskiid + alphaiid[ggl$y==1,]
	###
	var2 <- t(iidfull) %*% iidfull
	se <- cbind(diag(var2)^.5); colnames(se) <- "se"
} else { iidfull <- NULL; var2 <- NULL; se <- NULL} ## }}} 

var.naive=udca$robvar.gamma
se.naive=matrix(diag(var.naive)^.5,nrow(var.naive),1); 
colnames(se.naive) <- "se.naive"

res <- list(iid=iidfull,coef=udca$gamma,var.naive=var.naive,
	    se.naive=se.naive,var=var2,se=se,
	    comprisk.ipw=udca)
class(res) <- "comprisk.ipw"
return(res)
## }}} 
} 


### glmformula must have cause specific covariates if 
###      logit(P(e_i==1| V_i , cause_i==1))
###      logit(P(e_i==1| V_i , cause_i!=1))
###   potentially 
###      logit(P(e_i==1| V_i , cause_i==i)) can be specified where V_i are the covariates that 
###   are always observed and used for estimating the probability of sampling 
###   then glmformula must allow this but can still use comprisk.ipw function 


summary.comprisk.ipw <- function(object,digits=3,...)
{ ## {{{ 

tval <- object$coef/object$se
pval <- 2*(1-pnorm(abs(tval)))
res <- cbind(object$coef,object$se,object$se.naive,pval)
colnames(res) <- c("coef","se","se.naive","pval")

return(res)
} ## }}} 

coef.comprisk.ipw<- function(object,digits=3,...)
{ ## {{{ 
summary.comprisk.ipw(object)
} ## }}} 

print.comprisk.ipw  <-  function(x,...)
{  ## {{{ 
summary.comprisk.ipw(x)
}  ## }}} 


