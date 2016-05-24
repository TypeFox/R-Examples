feerror<-function(env, beta,sige, effects, method,rho, legacy){

	y<-get("y", envir = env)
	x<-get("x", envir = env)
	

	yt<-get("yt", envir = env)
	xt<-get("xt", envir = env)	
	N<-get("n", envir = env)
	T<-get("T", envir = env)
	NT<-get("NT", envir = env)
	k<-get("k", envir = env)
	listw<-get("listw", envir = env)
	inde<-get("inde", envir = env)


mx<-apply(x,2,mean)
intercept<-mean(y)-mx%*%beta
if (effects =="spfe"){
	ysms<-get("ysms", envir = env)
	xsms<-get("xsms", envir = env)

	sige<-as.numeric(sige)
	res.sfe <- as.matrix(ysms) - xsms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.sfe,T) + as.numeric(intercept)
	res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe<-sqrt(diag(res.var.sfe))
	res.t.sfe <- res.sfe / res.se.sfe 
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + N
	res.e <- y - xhat
	FE.out<-list(res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
	}
	
if (effects== "tpfe")	{
	
	ytms<-get("ytms", envir = env)
	xtms<-get("xtms", envir = env)

 	sige<-as.numeric(sige)
	res.tfe <- as.matrix(ytms) - xtms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.tfe,each=N) + as.numeric(intercept)
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe <-sqrt(diag(res.var.tfe))
	res.t.tfe <- res.tfe/res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + T
		res.e <- y - xhat
		FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
		
if (effects== "sptpfe"){
	
	ysms<-get("ysms", envir = env)
	xsms<-get("xsms", envir = env)

	ytms<-get("ytms", envir = env)
	xtms<-get("xtms", envir = env)

	sige<-as.numeric(sige)
	res.sfe <- as.matrix(ysms) - xsms %*% as.matrix(beta) - as.numeric(intercept)
	res.tfe <- as.matrix(ytms) - xtms %*% as.matrix(beta) - as.numeric(intercept)
		res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe <-sqrt(diag(res.var.sfe))
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe<-sqrt(diag(res.var.tfe))
	res.t.sfe <- res.sfe / res.se.sfe
	res.t.tfe <- res.tfe / res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	xhat<- x %*% as.matrix(beta) + rep(res.sfe,T) + rep(res.tfe,each=N) + as.numeric(intercept)
	N.vars <- k + N + T - 1
	res.e <- y - xhat
FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
		
if (effects=="pooled") {
	xhat <-   x %*% as.matrix(beta)
	res.e <- y - xhat
	FE.out<-list(xhat=xhat,N.vars=k,res.e=res.e)
	}

	yhat <- xhat
	ywhat <-  xt %*% beta
	r1 <- as.matrix(yt - mean(yt))
	r2 <- as.matrix(ywhat - mean(ywhat))
	r1r2 <- crossprod(r1,r2)
	r1r1 <- crossprod(r1)
	r2r2 <- crossprod(r2)
	res.corr <- as.numeric(r1r2^2) / (as.numeric(r1r1)*as.numeric(r2r2))
	
FE.out <- list(FE.out, res.corr=res.corr)
	}
	
#felag<-function(y,x,wy,ysms,xsms,ytms, xtms, wytms, wysms, beta,sige,yt,xt,N,T,NT,k,effects,method, rho,listw,inde){

felag<-function(env, beta,sige, effects, method, lambda, legacy, zero.policy){
	y<-get("y", envir = env)
	x<-get("x", envir = env)
	wy<-get("wy", envir = env)
	

	yt<-get("yt", envir = env)
	xt<-get("xt", envir = env)	
	N<-get("n", envir = env)
	T<-get("T", envir = env)
	NT<-get("NT", envir = env)
	k<-get("k", envir = env)
	listw<-get("listw", envir = env)
	inde<-get("inde", envir = env)
		
		mx<-apply(x,2,mean)
		intercept <- mean(y)- mean(wy)*lambda -  mx%*%beta

if (effects=="spfe"){
	ysms<-get("ysms", envir = env)
	xsms<-get("xsms", envir = env)
	wysms<-get("wysms", envir = env)

	res.sfe <- as.matrix(ysms) - as.matrix(wysms) *lambda - xsms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.sfe,T) + as.numeric(intercept)
	res.var.sfe<- (sige / T)  + diag((as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms))))
	res.se.sfe<-sqrt(res.var.sfe)
	res.t.sfe <- res.sfe / res.se.sfe 
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + N
	res.e <- y - xhat - lambda* wy
FE.out<-list(res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, 	res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
	}
	
if (effects== "tpfe")	{
	ytms<-get("ytms", envir = env)
	xtms<-get("xtms", envir = env)
	wytms<-get("wytms", envir = env)	

	res.tfe <- as.matrix(ytms) - as.matrix(wytms)* lambda - xtms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.tfe,each=N) + as.numeric(intercept)
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe <-sqrt(diag(res.var.tfe))
	res.t.tfe <- res.tfe/res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + T
	res.e <- y - xhat - lambda * wy
FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, intercept=intercept, 	res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects== "sptpfe"){
	ysms<-get("ysms", envir = env)
	xsms<-get("xsms", envir = env)
	wysms<-get("wysms", envir = env)

	ytms<-get("ytms", envir = env)
	xtms<-get("xtms", envir = env)
	wytms<-get("wytms", envir = env)	

	res.sfe <- as.matrix(ysms) - as.matrix(wysms) * lambda - xsms %*% as.matrix(beta) - as.numeric(intercept)
	res.tfe <- as.matrix(ytms) - as.matrix(wytms) * lambda - xtms %*% as.matrix(beta) - as.numeric(intercept)
	res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe <-sqrt(diag(res.var.sfe))
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe<-sqrt(diag(res.var.tfe))
	res.t.sfe <- res.sfe / res.se.sfe
	res.t.tfe <- res.tfe / res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	xhat<- x %*% as.matrix(beta) + rep(res.sfe,T) + rep(res.tfe,each=N) + as.numeric(intercept)
	N.vars <- k + N + T - 1
	res.e <- y - xhat - lambda * wy
FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects=="pooled") {
	xhat <-   x %*% as.matrix(beta)
	res.e <- y - xhat - lambda* wy
	FE.out<-list(xhat=xhat,N.vars=k,res.e=res.e)
	}

if(legacy){
	W <- listw2dgCMatrix(listw, zero.policy = zero.policy)
	IrW<- sparseMatrix(i=1:N, j=1:N, x=1) -lambda*W
	IrWi<-solve(IrW)
	xtb <- xt %*% beta
	yhat <- unlist(tapply(xhat,inde, function(u) IrWi %*% u))
	ywhat <- unlist(tapply(xtb,inde, function(u) IrWi %*% u))
	r1 <- as.matrix(yt - mean(yt))
	r2 <- as.matrix(ywhat - mean(ywhat))
	r1r2 <- crossprod(r1,r2)
	r1r1 <- crossprod(r1)
	r2r2 <- crossprod(r2)
	res.corr <- as.numeric(r1r2^2) / (as.numeric(r1r1)*as.numeric(r2r2))
}
else res.corr <- NULL  

FE.out <- list(FE.out, res.corr=res.corr)
FE.out
	}
	
	
	
effects.splm<-function(object,...){
	x<-object
	if (class(x) != "splm") stop(paste("methos not implemented for objects of class", class(x)))
	if (class(x) != "splm" && (x$type != "fixed effects lag" || x$type != "fixed effects error")) stop(paste("methos not implemented for objects of type", x$type))
	all.FE<-x$res.eff[[1]]
	effects <- x$effects
if (effects=="pooled") stop("No fixed effects available if effects == pooled")
if(effects=="spfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.se.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	SP.EFF <- all.FE$res.sfe
	se.SP.EFF <- all.FE$res.se.sfe 
	z <- SP.EFF/se.SP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   SETable <- cbind(SP.EFF,se.SP.EFF,z,p)
	colnames(SETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,SETable=SETable, effects=effects)
	}
if(effects=="tpfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.se.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	TP.EFF <- all.FE$res.tfe
	se.TP.EFF <- all.FE$res.se.tfe 
	z <- TP.EFF/se.TP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   TETable <- cbind(TP.EFF,se.TP.EFF,z,p)
	colnames(TETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,TETable=TETable,effects=effects)
	}
if(effects=="sptpfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.se.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	SP.EFF <- all.FE$res.sfe
	se.SP.EFF <- all.FE$res.se.sfe 
	z <- SP.EFF/se.SP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   SETable <- cbind(SP.EFF,se.SP.EFF,z,p)
	colnames(SETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
		TP.EFF <- all.FE$res.tfe
	se.TP.EFF <- all.FE$res.se.tfe 
	z <- TP.EFF/se.TP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   TETable <- cbind(TP.EFF,se.TP.EFF,z,p)
	colnames(TETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,SETable=SETable,TETable=TETable,effects=effects)
	}
res
class(res) <- "effects.splm"
return(res)
	}
	
	
print.effects.splm <-
function(x,digits= max(3, getOption("digits") - 2),
...){
	object<-x
	effects<-object$effects
if(effects=="tpfe"){
	  cat("\nIntercept:\n")
  printCoefmat(object$INTTable,digits=digits, signif.legend=FALSE)

 cat("\n")  

	  cat("\nTime period fixed effects:\n")
  printCoefmat(object$TETable,digits=digits)

out<-rbind(object$INTTable,object$TETable)
}
	
if(effects=="spfe"){
	  cat("\nIntercept:\n")
  printCoefmat(object$INTTable,digits=digits,signif.legend=FALSE)

 cat("\n")  
 
	  cat("\nSpatial fixed effects:\n")
  printCoefmat(object$SETable,digits=digits)


out<-rbind(object$INTTable,object$SETable)
}
if(effects=="sptpfe"){
	  cat("\nIntercept:\n")
  printCoefmat(object$INTTable,digits=digits,signif.legend=FALSE)

 cat("\n")  
 
	  cat("\nSpatial fixed effects:\n")
  printCoefmat(object$SETable,digits=digits,signif.legend=FALSE)
 
  cat("\n")  
   
  	  cat("\nTime period fixed effects:\n")
  printCoefmat(object$TETable,digits=digits)

out<-rbind(object$INTTable,object$SETable,object$TETable)
}
	}

