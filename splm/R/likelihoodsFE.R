#### SAR PANEL
`conclikpan` <- function(lambda, env){
	
	e0e0 <- get("e0e0", envir = env)
	e1e1 <- get("e1e1", envir = env)
	e0e1 <- get("e0e1", envir = env)
	NT <- get("NT", envir = env)
	T <- get("T", envir = env)
	
	Nsig <- e0e0 - 2*lambda*e0e1 + lambda*lambda*e1e1
	sigma2 <- Nsig/NT
	
	ldet <-  do_ldet(lambda, env)

ret <- - (NT/2)*log(Nsig)  + T * ldet  


	  if (get("verbose", envir=env)) 
        cat("lambda:\t", lambda, "\tfunction value:\t", ret, 
            "\n")
	ret
}


splaglm<-function(env, zero.policy = zero.policy, interval = interval, con = con, llprof = llprof, tol.solve= tol.solve, Hess = Hess, method = method, LeeYu = LeeYu, effects = effects){

xt <- get("xt", envir = env)
yt <- get("yt", envir = env)
wyt <- get("wyt", envir = env)
con<-get("con", envir = env)
NT<-get("NT", envir = env)
T<-get("T", envir = env)
n <- NT/T
listw<-get("listw", envir = env)
inde<-get("inde", envir = env)
interval1 <- get("interval1", envir = env)

        XpX<-crossprod(xt)
		b0<-solve(XpX,crossprod(xt,yt)) ####y on X
		b1<-solve(XpX,crossprod(xt,wyt)) ####Wy on x
		e0<-yt - xt%*% b0
		e1<-wyt - xt%*% b1
		e0e0<-crossprod(e0)
		e1e1<-crossprod(e1)
		e0e1<-t(e1)%*%e0


assign("e0e0", e0e0, envir = env)		
assign("e1e1", e1e1, envir = env)		
assign("e0e1", e0e1, envir = env)		
		
 
opt <- optimize(conclikpan,  interval = interval1, maximum = TRUE, env = env, tol = con$tol.opt)
#opt <- nlminb(0.02138744, conclikpan,  lower = interval[1], upper= interval[2],  env = env)

        lambda <- opt$maximum
		
    if (isTRUE(all.equal(lambda, interval[1])) || isTRUE(all.equal(lambda,interval[2]))) 
        warning("lambda on interval bound - results should not be used")

        names(lambda) <- "lambda"
        LL <- opt$objective
        optres <- opt

	lm.lag <- lm((yt - lambda * wyt) ~ xt - 1)
	p <- lm.lag$rank
    r <- residuals(lm.lag)
    fit <- yt - r
    names(r) <- names(fit)
    SSE <- crossprod(residuals(lm.lag))
    s2 <- as.numeric(SSE)/NT

if(LeeYu && effects == "spfe") s2 <- (T/(T-1)) * as.numeric(s2)	
if(LeeYu && effects == "tpfe") s2 <- (n/(n-1)) * as.numeric(s2)	

	betas <- coefficients(lm.lag)
	 # betas <- b0 - lambda*b1
	names(betas) <- colnames(xt)
	coefs <- c(lambda, betas)


if(LeeYu && effects == "sptpfe"){
	   
	    tr <- function(A) sum(diag(A))
        W <-listw2dgCMatrix(listw, zero.policy = zero.policy)
        A <- solve(sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1) - lambda * W)
        WA <- W %*% A
		lag <- function(q) trash<-unlist(tapply(q,inde,function(TT) as.matrix(WA %*% TT), simplify=TRUE))	  
		lag2 <- function(q) trash<-unlist(tapply(q,inde,function(TT) as.matrix(t(WA)%*%TT), simplify=TRUE))
		WAxt <- apply(as.matrix(xt),2,lag)
        WAWAxt<-apply(WAxt,2,lag2)
        xtWAWAxt <- crossprod(xt,WAWAxt)
        xtWAxt <- crossprod(xt,WAxt)
        xtxt <- crossprod(xt) 
	    one  <- T*(tr(WA %*% WA) + tr(t(WA) %*% WA))	    
        two <- 1/as.numeric(s2) * t(betas) %*% xtWAWAxt  %*% betas        
		V <- one + two
		zero <- rbind(rep(0, length(betas)))
        col1 <- rbind(NT/(2 * (s2^2)), T*tr(WA)/s2, t(zero))
        three <- (1/as.numeric(s2)) * xtWAxt %*% betas
        col2 <- rbind(T*tr(WA)/s2, V, three )
        col3 <- rbind(zero, t(three), 1/as.numeric(s2)* xtxt)
        asyvar <- cbind(col1, col2, col3)
        asyv <- solve(asyvar, tol = con$tol.solve)
		rownames(asyv) <- colnames(asyv) <- c("sigma","lambda", colnames(xt))

		init <- c((T/(T+1)), rep(1,p+1))	

		a3 <- rep(0,p)
		a2 <- 1/(1 - lambda)
		a1 <- 1/(2*s2)
		a <- c(a1,a2,a3)
		Bhat <- - (asyv/n) %*% a


		coefs1 <- c(s2,  lambda, betas)
		Theta2 <- init * coefs1 + Bhat
 		betas <- as.numeric(Theta2[(3:(2+p))])
 		names(betas) <- colnames(xt)  
 		lambda <-Theta2[2] 
        names(lambda) <- "lambda"
 		s2 <-  Theta2[1]
	    coefs <- c(lambda, betas)

	
}	

### add numerical hessian
if(Hess){
	
        fd <- fdHess(coefs, f_sarpanel_hess, env, LeeYu = LeeYu, effects = effects)
        mat <- fd$Hessian
		fdHess<- solve(-(mat), tol.solve = tol.solve)
        rownames(fdHess) <- colnames(fdHess) <- c("lambda", colnames(xt))
        
        lambda.se <- fdHess[1, 1]
        sig.se <- NULL
        asyvar1 <- vcov(lm.lag)
         rest.se<- NULL   
            
}

else{

        tr <- function(A) sum(diag(A))
        W <-listw2dgCMatrix(listw, zero.policy = zero.policy)
        A <- solve(sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1) - lambda * W)
        WA <- W %*% A
		lag <- function(q) trash<-unlist(tapply(q,inde,function(TT) as.matrix(WA %*% TT), simplify=TRUE))	  
		lag2 <- function(q) trash<-unlist(tapply(q,inde,function(TT) as.matrix(t(WA)%*%TT), simplify=TRUE))
		WAxt <- apply(as.matrix(xt),2,lag)
        WAWAxt<-apply(WAxt,2,lag2)
        xtWAWAxt <- crossprod(xt,WAWAxt)
        xtWAxt <- crossprod(xt,WAxt)
        xtxt <- crossprod(xt) 

if(LeeYu && effects == "spfe"){
	T <- T- 1
	NT <- n*T
}	

if(LeeYu && effects == "tpfe"){
	n <- n-1
	NT <- n*T
}	
        
 if(LeeYu && effects == "sptpfe"){
	n <- n-1
	T <- T-1
	NT <- n*T
}		
        one  <- T*(tr(WA %*% WA) + tr(t(WA) %*% WA))
        two <- 1/as.numeric(s2) * t(betas) %*% xtWAWAxt  %*% betas
		V <- one + two
		zero <- rbind(rep(0, length(betas)))
        col1 <- rbind(NT/(2 * (s2^2)), T*tr(WA)/s2, t(zero))
        three <- (1/as.numeric(s2)) * xtWAxt %*% betas
        col2 <- rbind(T*tr(WA)/s2, V, three )
        col3 <- rbind(zero, t(three), 1/as.numeric(s2)* xtxt)
        asyvar <- cbind(col1, col2, col3)
        asyva <- solve(asyvar, tol = con$tol.solve)
        rownames(asyva) <- colnames(asyva) <- c("sigma","lambda", colnames(xt))
        
        lambda.se <- asyva[2, 2]        
        rest.se <- sqrt(diag(asyva))[-c(1:2)]
        sig.se <- sqrt(asyva[1, 1]) 
        asyv <- asyva[-1,-1]      
        asyvar1 <- as.matrix(asyva[-c(1,2),-c(1,2)])
        rownames(asyvar1) <- colnames(asyvar1) <- colnames(xt)


}

if(Hess) asyv <- NULL        
else asyv <- asyv


 
    	return<-list(coeff = betas, lambda = lambda, s2 = s2, rest.se = rest.se, lambda.se = lambda.se, sig.se = sig.se, asyvar1 = asyvar1,  residuals = r, asyv = asyv)
} 



f_sarpanel_hess <- function (coefs, env, effects = effects, LeeYu = LeeYu) 
{
	
	T<-get("T", envir = env)
	NT<-get("NT", envir = env)
    n <- NT/T
    
if(LeeYu && effects == "spfe"){

	T <- T- 1
	NT <- n*T
}	

if(LeeYu && effects == "tpfe"){
	n <- n-1
	NT <- n*T
}	


if(LeeYu && effects == "sptpfe"){
	n <- n-1
	T <- T-1
	NT <- n*T
}		

    lambda <- coefs[1]
    beta <- coefs[-1]
    SSE <- sar_hess_sse_panel(lambda, beta, env)
    s2 <- SSE /n
    
    ldet <- do_ldet(lambda, env, which = 1)
    ret <- (T * ldet  - ((n*T/2) * log(2 * pi)) - (n*T/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir = env)) 
        cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet," SSE:", SSE, "\n")
    ret
}

sar_hess_sse_panel <- function (lambda, beta, env) 
{
    yl <- get("yt", envir = env) - lambda * get("wyt", envir = env) 
    res <- yl - (get("xt", envir = env) %*% beta)
    SSE <- c(crossprod(res))
    SSE
}


####### ERROR MODEL 
sarpanelerror<-function (rho, env=env) 
{
	yt<- get("yt", envir = env)
	xt<- get("xt", envir = env)
	wyt<- get("wyt", envir = env)
	wxt<- get("wxt", envir = env)
	wy<- get("wy", envir = env)
	wx<- get("wx", envir = env)

	listw<- get("listw", envir = env)
	NT<- get("NT", envir = env)
	inde<- get("inde", envir = env)
	T<- get("T", envir = env)
	
    yco <- yt - rho * wyt
    xco <- xt - rho * wxt
    bb<- solve(crossprod(xco),crossprod(xco, yco) )

    ehat<- yco - xco %*% bb
    SSE <- crossprod(ehat)
  ldet <- do_ldet(rho, env)

    ret <- T*ldet - (NT/2) * log(SSE) 

if (get("verbose", envir = env)) 
        cat("rho:", rho, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
 ret
}





sperrorlm <- function(env, zero.policy = zero.policy, interval = interval, con = con, llprof = llprof, tol.solve= tol.solve, Hess = Hess, LeeYu = LeeYu, effects = effects){

xt <- get("xt", envir = env)
yt <- get("yt", envir = env)
wyt <- get("wyt", envir = env)
wxt<-get("wxt", envir = env)

wy <- get("wy", envir = env)
wx<-get("wx", envir = env)

con<-get("con", envir = env)
NT<-get("NT", envir = env)
T<-get("T", envir = env)
n <- NT/T
listw<-get("listw", envir = env)
inde<-get("inde", envir = env)
interval <- get("interval1", envir = env)

opt <- optimize(sarpanelerror, interval = interval, maximum = TRUE, env = env, tol = con$tol.opt)


#opt <- nlminb(0.5,sarpanelerror,lower = interval[1], upper= interval[2], env = env)
#print(opt)

        rho <- opt$maximum
        names(rho) <- "rho"
        LL <- opt$objective

    if (isTRUE(all.equal(rho, interval[1])) || isTRUE(all.equal(rho,interval[2]))) 
        warning("rho on interval bound - results should not be used")

    lm.target <- lm(I(yt - rho * wyt) ~ I(xt - rho * wxt) - 
        1)
    r <- as.vector(residuals(lm.target))
    p <- lm.target$rank
    SSE <- crossprod(residuals(lm.target))
    s2 <- as.numeric(SSE)/NT

if(LeeYu && effects == "spfe") s2 <- (T/(T-1)) * as.numeric(s2)	
if(LeeYu && effects == "tpfe") s2 <- (n/(n-1)) * as.numeric(s2)	

    rest.se <- (summary(lm.target)$coefficients[, 2]) * sqrt((NT - p)/NT)     
    betas <- coefficients(lm.target)
    names(betas) <- colnames(xt)  
     coefs <- c(rho, betas) 

if(LeeYu && effects == "sptpfe"){
	    
	    tr <- function(A) sum(diag(A))
        W <- listw2dgCMatrix(listw, zero.policy = zero.policy)
        A <- solve(sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1)  - rho * W)
        WA <- W %*% A
        asyvar <- matrix(0, nrow = 2 + p, ncol = 2 + p)
        asyvar[1, 1] <- NT/(2 * (s2^2))
        asyvar[2, 1] <- asyvar[1, 2] <- T*tr(WA)/s2
        asyvar[2, 2] <- T*(tr(WA %*% WA) + tr(t(WA) %*% WA))
        asyvar[3:(p + 2), 3:(p + 2)] <- 1/as.numeric(s2) * (t(xt - rho *wxt) %*% (xt - rho * wxt)) 
        asyv <- solve(asyvar, tol = con$tol.solve)
        rownames(asyv) <- colnames(asyv) <- c("sigma","rho", colnames(xt))

        s2.se <- sqrt(asyv[1, 1])
        rho.se <- asyv[2, 2]
        asyvar1 <- asyv[-c(1,2),-c(1,2)]
		init <- c((T/(T+1)), rep(1,p+1))	

		a3 <- rep(0,p)
		a2 <- 1/(1 - rho)
		a1 <- 1/(2*s2)
		a <- c(a1,a2,a3)
		Bhat <- - (asyv/n) %*% a


		coefs1 <- c(s2, rho, betas)
		Theta2 <- init * coefs1 + Bhat
 		betas <- as.numeric(Theta2[(3:(2+p))])
 		names(betas) <- colnames(xt)  
 		rho <-Theta2[2] 
        names(rho) <- "rho"
 		s2 <-  Theta2[1]
	    coefs <- c(rho, betas)

	
	
	
}

if(Hess){
	
	    fd <- fdHess(coefs, sarpanelerror_hess, env, LeeYu = LeeYu, effects = effects)
        mat <- fd$Hessian
		fdHess<- solve(-(mat), tol.solve = tol.solve)
        rownames(fdHess) <- colnames(fdHess) <- c("rho",colnames(xt))

            rho.se <- fdHess[1, 1]
            s2.se <- NULL
            asyvar1 <- vcov(lm.target)
}
else{
    
        tr <- function(A) sum(diag(A))
        W <- listw2dgCMatrix(listw, zero.policy = zero.policy)
        A <- solve(sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1)  - rho * W)
        WA <- W %*% A

if(LeeYu && effects == "spfe"){

	T <- T- 1
	NT <- n*T

}	

if(LeeYu && effects == "tpfe"){

	n <- n-1
	NT <- n*T

}	

 if(LeeYu && effects == "sptpfe"){
	n <- n-1
	T <- T-1
	NT <- n*T
}		

        asyvar <- matrix(0, nrow = 2 + p, ncol = 2 + p)
        asyvar[1, 1] <- NT/(2 * (s2^2))
        asyvar[2, 1] <- asyvar[1, 2] <- T*tr(WA)/s2
        asyvar[2, 2] <- T*(tr(WA %*% WA) + tr(t(WA) %*% WA))
        asyvar[3:(p + 2), 3:(p + 2)] <- 1/as.numeric(s2) * (t(xt - rho *wxt) %*% (xt - rho * wxt)) 
        asyva <- solve(asyvar, tol = con$tol.solve)
        rownames(asyva) <- colnames(asyva) <- c("sigma","rho", colnames(xt))
        s2.se <- sqrt(asyva[1, 1])
        rho.se <- asyva[2, 2]
        asyvar1 <- asyva[-c(1,2),-c(1,2)]
        asyv <- asyva[-1,-1]
        rownames(asyvar1) <- colnames(asyvar1) <- colnames(xt)

        

}

if(Hess) asyv <- NULL        
else asyv <- asyv


	return<-list(coeff=betas, rho = rho, s2 = s2, rest.se = rest.se, rho.se = rho.se, s2.se = s2.se, asyvar1=asyvar1, residuals = r, asyv = asyv)
}



sarpanelerror_hess<-function (coef, env=env, LeeYu = LeeYu, effects = effects) 
{
	yt<- get("yt", envir = env)
	xt<- get("xt", envir = env)
	wyt<- get("wyt", envir = env)
	wxt<- get("wxt", envir = env)
	wy<- get("wy", envir = env)
	wx<- get("wx", envir = env)

	listw<- get("listw", envir = env)
	NT<- get("NT", envir = env)
	inde<- get("inde", envir = env)
	T<- get("T", envir = env)
	n <- NT/T

if(LeeYu && effects == "spfe"){
	T <- T- 1
	NT <- n*T
}	

if(LeeYu && effects == "tpfe"){
	n <- n-1
	NT <- n*T
}	

 if(LeeYu && effects == "sptpfe"){
	n <- n-1
	T <- T-1
	NT <- n*T
}		

	rho <- coef[1]
	bb <- coef[-1]
	 
     yco <- yt - rho * wyt
     xco <- xt - rho * wxt

     ehat<- yco - xco %*% bb
    SSE <- crossprod(ehat)

  ldet <- do_ldet(rho, env)

    ret <- T*ldet - (NT/2) * log(SSE) 

if (get("verbose", envir = env)) 
        cat("rho:", rho, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
 ret
}


###SARAR MODEL

sacsarpanel<-function (coefs, env){

	lambda <- coefs[1]
    rho <- coefs[2]
  	 T<-get("T", envir = env)
    n <- get("n", envir = env)

    SSE <- sacsarpanel_sse(coefs, env)
    s2 <- SSE/n
    ldet1 <- do_ldet(lambda, env, which = 1)
    ldet2 <- do_ldet(rho, env, which = 2)

ret <- (T * ldet1 + T * ldet2 - (((n*T)/2) * (log(2 * pi)+1)) - (n*T/2) * log(s2))
                        # - (1/(2 * (s2))) * SSE)
if(get("verbose", envir = env)) cat("lambda:", lambda, " rho:", rho, " function:", 
             ret, " Jacobian1:", ldet1, " Jacobian2:", ldet2, 
             " SSE:", SSE, "\n")
-ret
}


sacsarpanel_sse <- function (coefs, env) 
{
    lambda <- coefs[1]
    rho <- coefs[2]
    yl <- get("yt", envir = env) - lambda * get("wyt", envir = env) - 
        rho * get("w2yt", envir = env) + rho * lambda * get("w2wyt", 
        envir = env)
    xl <- get("xt", envir = env) - rho * get("wxt", envir = env)
    xl.q <- qr.Q(qr(xl, LAPACK = get("LAPACK", envir = env)))
    xl.q.yl <- crossprod(xl.q, yl)
    SSE <- crossprod(yl) - crossprod(xl.q.yl)
    SSE
}


spsararlm<-function(env, zero.policy = zero.policy, con = con, llprof = llprof, tol.solve= tol.solve, Hess = Hess, method = method, LeeYu = LeeYu, effects = effects){

xt <- get("xt", envir = env)
yt <- get("yt", envir = env)
wyt <- get("wyt", envir = env)
w2yt <- get("w2yt", envir = env)
w2wyt <- get("w2wyt", envir = env)
wxt<-get("wxt", envir = env)

wy <- get("wy", envir = env)
wx<-get("wx", envir = env)
	
NT<-get("NT", envir = env)
T<-get("T", envir = env)
n<-get("n", envir = env)

listw<-get("listw", envir = env)
listw2<-get("listw2", envir = env)
inde<-get("inde", envir = env)
interval1 <- get("interval1", envir = env)
interval2 <- get("interval2", envir = env)
	
    pars <- con$pars
    lower <- c(interval1[1], interval2[1])
    upper <- c(interval1[2], interval2[2])

    if (!is.null(llprof)) {
        llrho <- NULL
        lllambda <- NULL
        if (length(llprof) == 1L) {
            llrho <- seq(lower[2], upper[2], length.out = llprof)
            lllambda <- seq(lower[1], upper[1], length.out = llprof)
            llprof <- as.matrix(expand.grid(lllambda, llrho))
        }
        ll_prof <- numeric(nrow(llprof))
        for (i in 1:nrow(llprof)) ll_prof[i] <- sacsarpanel(llprof[i, 
            ], env = env)
        # nm <- paste(method, "profile", sep = "_")
        # timings[[nm]] <- proc.time() - .ptime_start
        # .ptime_start <- proc.time()
    }
    if (is.null(pars)) {
        if (con$npars == 4L) {
            xseq <- c(lower[1], 0, upper[1], upper[1]) * 0.8
            yseq <- c(upper[2], 0, upper[2], lower[2]) * 0.8
            mpars <- cbind(xseq, yseq)
        }
        else {
            xseq <- seq(lower[1], upper[1], (upper[1] - lower[1])/2) * 
                0.8
            yseq <- seq(lower[2], upper[2], (upper[2] - lower[2])/2) * 
                0.8
            mpars <- as.matrix(expand.grid(xseq, yseq))
        }
    }
    else {
        mxs <- NULL
    }
    if (con$opt_method == "nlminb") {
        if (is.null(pars)) {
            mxs <- apply(mpars, 1, function(pp) -nlminb(pp, sacsarpanel, 
                env = env, control = con$opt_control, lower = lower, 
                upper = upper)$objective)
            pars <- mpars[which.max(mxs), ]
            optres <- nlminb(pars, sacsarpanel, env = env, control = con$opt_control, 
                lower = lower, upper = upper)
        }
        else {
            optres <- nlminb(pars, sacsarpanel, env = env, control = con$opt_control, 
                lower = lower, upper = upper)
        }
    }
    else if (con$opt_method == "L-BFGS-B") {
        if (is.null(pars)) {
            mxs <- apply(mpars, 1, function(pp) -optim(pars, 
                sacsarpanel, env = env, method = "L-BFGS-B", control = con$opt_control, 
                lower = lower, upper = upper)$objective)
            pars <- mpars[which.max(mxs), ]
            optres <- optim(pars, sacsarpanel, env = env, method = "L-BFGS-B", 
                control = con$opt_control, lower = lower, upper = upper)
        }
        else {
            optres <- optim(pars, sacsarpanel, env = env, method = "L-BFGS-B", 
                control = con$opt_control, lower = lower, upper = upper)
        }
    }
    else {
        if (is.null(pars)) {
            mxs <- apply(mpars, 1, function(pp) -optim(pars, 
                sacsarpanel, env = env, method = con$opt_method, 
                control = con$opt_control)$objective)
            pars <- mpars[which.max(mxs), ]
            optres <- optim(pars, sacsarpanel, env = env, method = con$opt_method, 
                control = con$opt_control)
        }
        else {
            optres <- optim(pars, sacsarpanel, env = env, method = con$opt_method, 
                control = con$opt_control)
        }
    }
    
    LL <- -optres$objective
    if (optres$convergence != 0) 
        warning(paste("convergence failure:", optres$message))
	
	# print(optres)
    rho <- optres$par[2]
    names(rho) <- "rho"
    lambda <- optres$par[1]
    names(lambda) <- "lambda"


    
    lm.target <- lm(I(yt - lambda * wyt - rho * w2yt + rho * lambda * 
        w2wyt) ~ I(xt - rho * wxt) - 1)

    r <- as.vector(residuals(lm.target))
    fit <- as.vector(yt - r)
    p <- lm.target$rank
    SSE <- crossprod(residuals(lm.target))
    s2 <- as.numeric(SSE)/NT

if(LeeYu && effects == "spfe") s2 <- (T/(T-1)) * as.numeric(s2)	
if(LeeYu && effects == "tpfe") s2 <- (n/(n-1)) * as.numeric(s2)	

    betas <- coefficients(lm.target)
    names(betas) <- colnames(xt)  
	 coefs <- c(lambda, rho, betas)

if(LeeYu && effects == "sptpfe"){
	    
	    tr <- function(A) sum(diag(A))
        W1 <- listw2dgCMatrix(listw, zero.policy = zero.policy)
        W2 <- listw2dgCMatrix(listw2, zero.policy = zero.policy)
        Sl <- sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1)  - lambda * W1
        Rr <- sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1)  - rho * W2
        Slinv <- solve(Sl)
        Rrinv <- solve(Rr)
        Gmat <- W1 %*% Slinv
        Hmat <- W2 %*% Rrinv        
        It <- sparseMatrix(i=1:T, j=1:T, x=1)         
        Jn <- Diagonal(n) - (1/n) * outer(rep(1,n),rep(1,n))


# Equation 53 Lee and Yu         
        Wdot <- Rr %*% W1  %*% Rrinv
        Gdot <- Rr %*% Gmat  %*% Rrinv
        GS <- t(Gdot) + Gdot
        HS <- t(Hmat) + Hmat
        Rrbig <- kronecker(It,Rr)
        RriB  <- kronecker(It,Rrinv) 
        GdotB<-  kronecker(It,Gdot)
        WdotB<-  kronecker(It,Wdot)
        JnB <- kronecker(It,Jn)
        Xdot <-  Rrbig %*% xt  
        JXdot <- JnB %*% Xdot
        GdXdb <- GdotB %*% Xdot %*% betas
        JGdXdb <- JnB %*% GdotB %*% Xdot %*% betas
        
	    fp   <- (1/s2) * crossprod(GdXdb, JGdXdb)
        lala <- fp + (T * tr(GS %*% Jn %*% Gdot))
        laro <- T * tr(HS %*% Jn %*% Gdot) 
        lasi <- (1/s2) * T * tr(Gdot) 
        roro <- T * tr(HS %*% Hmat)
        rosi <- (1/s2) * T * tr(Hmat)
        sisi <- NT/(2*s2*s2)
        bebe <- (1/s2) * crossprod(Xdot, JXdot)       
        bela <- (1/s2) * crossprod(GdXdb, JXdot)
        
        asyvar <- matrix(0, nrow = 3 + p, ncol = 3 + p)
        asyvar[1:p, 1:p] <- as.matrix(bebe) 
        asyvar[p+1, 1:p] <- asyvar[1:p, p+1] <- as.numeric(bela)
        asyvar[p+2, 1:p] <- asyvar[1:p, p+2] <- 0
        asyvar[p+3, 1:p] <- asyvar[1:p, p+3] <- 0
        asyvar[p+2, p+1] <- asyvar[p+1, p+2] <- as.numeric(laro)
        asyvar[p+3, p+1] <- asyvar[p+1, p+3] <- as.numeric(lasi)
        asyvar[p+3, p+2] <- asyvar[p+2, p+3] <- as.numeric(rosi)        
        asyvar[1+p, 1+p] <- as.matrix(lala)
        asyvar[2+p, 2+p] <- as.matrix(roro)
        asyvar[3+p, 3+p] <- as.matrix(sisi)
        asyv <- solve(asyvar, tol = con$tol.solve)

		a1 <- rep(0,p)
		a2 <- as.numeric((1/n) * rep(1,n) %*% Gdot %*% rep(1,n))
		a3 <- as.numeric((1/n) * rep(1,n) %*% Hmat %*% rep(1,n))
		a4 <- as.numeric(1/(2*s2))
		a <- c(a1,a2,a3,a4)
		Bhat <- - asyv %*% a
		At <- matrix(0, nrow = 3 + p, ncol = 3 + p)
		At[(1:(p+2)), (1:(p+2))]<- diag((p+2))
		At[(1:(p+2)), 3+p] <- rep(0,p+2)
		At[3+p,(1:(p+2))] <- rep(0,p+2)
		At[3+p, 3+p] <- T/(T-1)
		coefs1 <- c(betas, lambda, rho, s2) 
		Theta1 <- coefs1 - (Bhat/n)
		Theta2 <- At %*% Theta1
 		betas <- Theta2[1:p]
 		names(betas) <- colnames(xt)  
 		lambda <-Theta2[p+1] 
 		rho <- Theta2[p+2]
 		names(rho) <- "rho"
        names(lambda) <- "lambda"
 		s2 <-  Theta2[p+3]
	    coefs <- c(lambda, rho, betas)
	
}


###Add the vc matrix exact
if(Hess){    

# if(LeeYu && effects == "sptpfe") stop("Numerical Hessian should not be calculated when 'LeeYu = TRUE' and effects are 'twoways' ")
	    
        fd <- fdHess(coefs, f_sacpanel_hess, env, LeeYu = LeeYu, effects = effects)
        mat <- fd$Hessian
		  fdHess<- solve(-(mat), tol.solve = tol.solve)
        rownames(fdHess) <- colnames(fdHess) <- c("lambda", "rho",colnames(xt))
            
            rho.se <- fdHess[2,2]
            lambda.se <- fdHess[1,1]
            asyvar1 <- vcov(lm.target)
            s2.se <- NULL
            }
            
            else{
            	
   
        tr <- function(A) sum(diag(A))
        W1 <- listw2dgCMatrix(listw, zero.policy = zero.policy)
        W2 <- listw2dgCMatrix(listw2, zero.policy = zero.policy)
        Sl <- sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1)  - lambda * W1
        Rr <- sparseMatrix(i=1:(NT/T), j=1:(NT/T), x=1)  - rho * W2
        Slinv <- solve(Sl)
        Rrinv <- solve(Rr)
        Gmat <- W1 %*% Slinv
        Hmat <- W2 %*% Rrinv        
        It <- sparseMatrix(i=1:T, j=1:T, x=1)         
        

# Equation 39 Lee and Yu         
        Wdot <- Rr %*% W1  %*% Rrinv
        Gdot <- Rr %*% Gmat  %*% Rrinv
        GS <- t(Gdot) + Gdot
        HS <- t(Hmat) + Hmat
        Rrbig <- kronecker(It,Rr)
        RriB  <- kronecker(It,Rrinv) 
        GdotB<-  kronecker(It,Gdot)
        WdotB<-  kronecker(It,Wdot)
        Xdot <-  Rrbig %*% xt  
        GdXdb <- GdotB %*% Xdot %*% betas
   
        
if(LeeYu && effects == "spfe"){
	T <- T- 1
	NT <- n*T
}	

if(LeeYu && effects == "tpfe"){
	n <- n-1
	NT <- n*T
}	

if(LeeYu && effects == "sptpfe"){
	n <- n-1
	T <- T-1
	NT <- n*T
}		
        fp   <- (1/s2) *crossprod(GdXdb)
        lala <- fp + (T * tr(GS %*% Gdot))
        laro <- T * tr(HS %*% Gdot) 
        lasi <- (1/s2) * T * tr(Gdot) 
        roro <- T * tr(HS %*% Hmat)
        rosi <- (1/s2) * T * tr(Hmat)
        sisi <- NT/(2*s2*s2)
        bebe <- (1/s2) * crossprod(Xdot)       
        bela <- (1/s2) * crossprod(GdXdb, Xdot)
        
        asyvar <- matrix(0, nrow = 3 + p, ncol = 3 + p)
        asyvar[1:p, 1:p] <- as.matrix(bebe) 
        asyvar[p+1, 1:p] <- asyvar[1:p, p+1] <- as.numeric(bela)
        asyvar[p+2, 1:p] <- asyvar[1:p, p+2] <- 0
        asyvar[p+3, 1:p] <- asyvar[1:p, p+3] <- 0
        asyvar[p+2, p+1] <- asyvar[p+1, p+2] <- as.numeric(laro)
        asyvar[p+3, p+1] <- asyvar[p+1, p+3] <- as.numeric(lasi)
        asyvar[p+3, p+2] <- asyvar[p+2, p+3] <- as.numeric(rosi)        
        asyvar[1+p, 1+p] <- as.matrix(lala)
        asyvar[2+p, 2+p] <- as.matrix(roro)
        asyvar[3+p, 3+p] <- as.matrix(sisi)
        asyva <- solve(asyvar, tol = con$tol.solve)
        rownames(asyva) <- colnames(asyva) <- c(colnames(xt), "lambda", "rho", "sigma")

        s2.se <- asyva[3+p, 3+p]
        rho.se <- asyva[2+p, 2+p]
        lambda.se <- asyva[1+p, 1+p]
        rest.se <- sqrt(diag(asyva))[-((p+1):(p+3))]
        asyvar1 <- asyva[-((p+1):(p+3)),-((p+1):(p+3))]
        asyv <- asyva[-(p+3),-(p+3)]


            	}

if(Hess) asyv <- NULL        
else asyv <- asyv
        

return<-list(coeff = betas, lambda = lambda, rho = rho, s2 = s2, asyvar1 = asyvar1, lambda.se = lambda.se, rho.se = rho.se, s2.se = s2.se, residuals = r, asyv = asyv)	
	}

f_sacpanel_hess <- function (coefs, env, LeeYu = LeeYu, effects = effects) 
{
	T<-get("T", envir = env)
	NT<-get("NT", envir = env)
	n<-get("n", envir = env)

if(LeeYu && effects == "spfe"){
	T <- T- 1
	NT <- n*T
}	

if(LeeYu && effects == "tpfe"){
	n <- n-1
	NT <- n*T
}	

if(LeeYu && effects == "sptpfe"){
	n <- n-1
	T <- T-1
	NT <- n*T
}		

    lambda <- coefs[1] 
    rho <- coefs[2]
    beta <- coefs[-(1:2)]
      SSE <- sar_sac_hess_sse_panel(lambda, rho, beta, env)
    # SSE <- sar_sac_hess_sse_panel(lambda, rho, beta, env)
    n <- NT/T
    # SSE<- s2 *n
     s2<- SSE / n
    ldet1 <- do_ldet(lambda, env, which = 1)
    ldet2 <- do_ldet(rho, env, which = 2)
   
#ret <- (T * ldet1 + T * ldet2 - (((n*T)/2) * (log(2 * pi))) - (n*T/2) * log(s2))
                        # - (1/(2 * (s2))) * SSE)
ret <- (T * ldet1 + T * ldet2 - ((n*T/2) * log(2 * pi)) - (n*T/2) * log(s2) - 
        (1/(2 * s2)) * SSE)


    if (get("verbose", envir = env)) cat("rho:", rho, "lambda:", lambda, " function:", ret, 
            " Jacobian1:", ldet1, " Jacobian2:", ldet2, " SSE:", 
            SSE, "\n")
    ret
}

sar_sac_hess_sse_panel <- function (lambda, rho,  beta, env) 
{
    yl <- get("yt", envir = env) - lambda * get("wyt", envir = env) - 
        rho * get("w2yt", envir = env) + rho * lambda * get("w2wyt", 
         envir = env)
         
    xl <- get("xt", envir = env) - rho * get("wxt", envir = env)
    res <- yl - (xl %*% beta)
    SSE <- c(crossprod(res))
    SSE
}




