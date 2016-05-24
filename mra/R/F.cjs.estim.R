F.cjs.estim <- function(capture, survival, histories, cap.init, sur.init, 
    group, nhat.v.meth=1, 
    c.hat=-1.0, df=NA, intervals=rep(1,ncol(histories)-1), 
    conf=0.95, link="logit", 
    control=mra.control() ){

start.sec <- proc.time()[3]

run.date <- Sys.time()

# ---- Verify some of the inputs
if( missing(histories) ){
    stop("Capture histories matrix must be specified")
}
if( missing(capture) ){
    stop("Capture covariates must be specified")
}
if( missing(survival) ){
    stop("Survival covariates must be specified")
}

if( length(union( unique(histories), c(0,1,2))) > 3 ) stop("Capture histories must consist of 0's, 1's, and 2's only.")

if( !is.numeric(nhat.v.meth) || !(nhat.v.meth %in% c(1,2,3)) ){
 	stop("value of 'nhat.v.meth' must be 1, 2, or 3")
}

# Validity of cov.meth and algorithm and other inputs is determined in mra.control()

# ---- Initialize some variables
hist.name <- deparse(substitute(histories))
cr.call <- match.call()

nan <- nrow( histories )
ns <- ncol( histories )

if( length(intervals) < (ns-1)){
    stop(paste("Too few time intervals specified. INTERVALS vector should have length", ns-1))
} else if(length(intervals) >= (ns-1)){
    intervals <- intervals[1:(ns-1)]
    intervals <- c(intervals, 0)  # Make this vector length ns, but we never use the last element in estimation.
}


# ---- Get the X and Y matricies.  After this, the x and Y matricies are in huge NAN by (ns*nx) matricies.
covars <- F.cr.model.matrix( capture, survival, nan, ns )  

nx <- covars$n.cap.covars  #nx and ny include intercept.  This is total number of parameters
ny <- covars$n.sur.covars


#   Set initial values if missing or short
if( missing(cap.init) ){
    cap.init <- rep(0,nx)
} else if(length(cap.init) < (nx) ){
    cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 

if( missing(sur.init) ){
    sur.init <- rep(0,ny)
} else if(length(sur.init) < (ny) ){
    sur.init <- c(sur.init, rep(0, ny-length(sur.init)))
} 

#   Set up the tolerance vector, if not specified, or if not long enough
if( length(control$tol) < (nx+ny) ){
    control$tol <- rep(control$tol, trunc((nx+ny) / length(control$tol))+1)[1:(nx+ny)]
} else if( length(control$tol > (nx+ny)) ){
    control$tol <- control$tol[1:(nx+ny)]
}


#   Fix up the group variable
if( missing( group )){
    group <- rep(1, nan)
    ng <- 1
} else {
    ng <- length(unique(group))
}

#   Transfer over c.hat
vif <- c.hat

#   Do the estimation, but first allocate room for answers
loglik <- deviance <- aic <- qaic <- chisq.vif <- df.vif <- 0
parameters <- se.param <- rep(0, nx + ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- s.hat <- se.s.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- rep(0, ns)
# on entry to .Fortran, maxfn is maximum number of function evals.  On return, maxfn is actual number of evaluations.
maxfn <- control$maxfn  
if( is.na(df) ){
    df.estimated <- 1  # Have MRAWIN estimate rank of var-covar matrix
} else {
    #df.estimated <- 0  # Don't bother, df either set by user or will use nx+ny
    df.estimated <- 1  #  Have MRAWIN estimate number of parameters.  Work out the one the user wants later
}


#   Re-code the link specification to integers
if( link=="logit" ){
    link.code <- 1
} else if( link == "sine" ){
    link.code <- 2
} else if( link == "hazard" ){
    link.code <- 3
} else {
    stop("Unknown link function specified.")
}

if(control$trace) cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")


ans <- .Fortran( "cjsmod", 
        nan         = as.integer(nan), 
        ns          = as.integer(ns), 
        nx          = as.integer(nx), 
        ny          = as.integer(ny), 
        ng          = as.integer(ng), 
        histories   = as.integer(histories), 
        group       = as.integer(group), 
        algorithm   = as.integer(control$algorithm), 
        cov.meth    = as.integer(control$cov.meth), 
        trace       = as.integer(control$trace),
        link        = as.integer(link.code),
        nhat.v.meth = as.integer(nhat.v.meth), 
        capX        = as.double(covars$capX), 
        survX       = as.double(covars$survX), 
        cap.init    = as.double(cap.init), 
        sur.init    = as.double(sur.init), 
        maxfn       = as.integer(maxfn),
        beta.tol.vec= as.double(control$tol), 
        loglik      = as.double(loglik), 
        deviance    = as.double(deviance), 
        aic         = as.double(aic), 
        qaic        = as.double(qaic), 
        vif         = as.double(vif), 
        chisq.vif   = as.double(chisq.vif), 
        df.vif      = as.double(df.vif), 
        parameters  = as.double(parameters),
        se.param    = as.double(se.param), 
        covariance  = as.double(covariance), 
        p.hat       = as.double(p.hat), 
        se.p.hat    = as.double(se.p.hat), 
        s.hat       = as.double(s.hat), 
        se.s.hat    = as.double(se.s.hat), 
        n.hat       = as.double(n.hat), 
        se.n.hat    = as.double(se.n.hat), 
        exit.code   = as.integer(exit.code), 
        cov.code    = as.integer(cov.code), 
        df.estimated= as.integer(df.estimated), 
        intervals   = as.double(intervals), 
        PACKAGE="mra" 
        ) 


if( control$trace ) cat(paste("Returned from MRAWIN. Fitting details written to MRA.LOG.\n", sep=""))


## ----- Reset missing standard errors to NA
ans$se.param[ ans$se.param < 0 ] <- NA

# ----- R does not preserve the matrix structures in .Fortran call.  Put matricies, 
#   which are now vectors, back to matricies.
covariance <- matrix( ans$covariance, nrow=nx+ny ) 
ans$p.hat      <- matrix( ans$p.hat, nrow=nan )
ans$se.p.hat   <- matrix( ans$se.p.hat, nrow=nan )
ans$s.hat      <- matrix( ans$s.hat, nrow=nan )
ans$se.s.hat   <- matrix( ans$se.s.hat, nrow=nan )

# ----- Work out exit codes
if( ans$exit.code==0 ){
    exit.mess = "FAILURE: Initial Hessian not positive definite"
} else if( ans$exit.code == 1 ){
    exit.mess = "SUCCESS: Convergence criterion met"
} else if( ans$exit.code == 2 ){
    exit.mess = "FAILURE: G'dX > 0, rounding error"
} else if( ans$exit.code == 3 ){
    exit.mess = "FAILURE: Likelihood evaluated too many times"
} else if( ans$exit.code == -1 ){
    exit.mess = "FAILURE: Algorithm 2 not implimented yet.  Contact Trent McDonald."
} else {
    exit.mess = "Unknown exit code"
}

if(control$algorithm == 1){
    alg.mess <- "Optimization by VA09AD."
} else {
    alg.mess <- "Unknown optimization routine."
}

if(control$cov.meth == 1){
    cov.mess = "Covariance from numeric derivatives."
} else if (control$cov.meth == 2){
    cov.mess = "Covariance from optimization Hessian."
} else {
    cov.mess = "Unkown covariance method."
}

if(ans$cov.code == 0){
    cov.mess = paste( cov.mess,  "SUCCESS: Non-singular covariance matrix.")
} else if( ans$cov.code == 1) {
    cov.mess = paste( cov.mess,  "ERROR: COVARIANCE MATRIX IS SINGULAR.")
} else {
    cov.mess = paste( cov.mess,  "ERROR: NON-NEGATIVE DEFINITE COVARIANCE MATRIX.")
}

# do some printouts
if( ans$exit.code != 1 ){
	cat(exit.mess)
	cat("\n")
} else if( ans$cov.code != 0 ){
	cat(cov.mess)
	cat("\n")
} else if( control$trace ){
	cat(paste( alg.mess, "\n", sep="" ))
	cat(paste( exit.mess, "\n", sep="" ))
	cat(paste( cov.mess, "\n", sep="" ))
}

# ----- Fix up capture and survival estimates
#   Wipe out the first capture probability and the last survival probability. These are computed 
#   by mrawin5 if there are covariates out there, but they are not used in the likelihood, and 
#   we can't really estimate them.
ans$p.hat[,1] <- NA
ans$s.hat[,ns] <- NA
ans$se.s.hat[,ns] <- NA
ans$se.p.hat[,1] <- NA

capcoef <- ans$parameters[1:nx]
se.capcoef <- ans$se.param[1:nx]
surcoef <- ans$parameters[(nx+1):(nx+ny)]
se.surcoef <- ans$se.param[(nx+1):(nx+ny)]

names(capcoef) <- covars$cap.vars
names(se.capcoef) <- covars$cap.vars
names(surcoef) <- covars$sur.vars
names(se.surcoef) <- covars$sur.vars


dimnames(covariance) <- list( c( paste("cap.",names( capcoef ),sep=""), paste("sur.",names( surcoef ), sep="")),
    c( paste("cap.",names( capcoef ),sep=""), paste("sur.",names( surcoef ), sep="")))

# ----- Fix up number of parameters.
if( is.na(df) ){
    df <- ans$df.estimated   # use the rank estimated by MRAWIN
} else if( df <= 0 ){
    df <- nx+ny  # assume full rank
} # else use the values supplied by user (unchanged from input)

# ----- Now that df is computed, recompute fit statistics
aic <- -2*ans$loglik + 2*df
qaic <- -2*(ans$loglik/ans$vif) + 2*df
aicc <- aic + (2*df*(df+1))/(nan - df - 1)
qaicc <- qaic + (2*df*(df+1))/(nan - df - 1)


# ----- LEAVE THIS CODE HERE FOR FUTURE REFERENCE 
# ----- Compute EBC of Peterson (1986), Stat & Prob letters, p227
#qf <- t(parameters) %*% covariance %*% parameters
#q1 <- qf/((df)*nan)
#q2 <- qf/(df)
#ebc <- -2*loglik + (df)*log(nan) + 
#    (df)*log( max(c( 1/nan, q1 ))) +
#    (df)*min( c( 1, q2 ) )
    
# ----- Put all the "auxillary" info together
aux <- list( call=cr.call, 
        nan=nan, 
        ns=ns, 
        nx=nx, 
        ny=ny, 
        link=link,
        cov.name=c(names(capcoef), names(surcoef)), 
        ic.name=hist.name, 
        mra.version=packageDescription("mra")$Version, 
        R.version=R.version.string,
        run.date=run.date )


# ----- Fix up the estimates of N.
names(ans$n.hat) <- dimnames(histories)[[2]]
num.observed <- c( t( rep(1, nrow(histories)) ) %*% (histories >= 1) )
crit.value <- qnorm( 1-((1-conf)/2) )
lower.ci <- ans$n.hat - crit.value * ans$se.n.hat
lower.ci <- ifelse(lower.ci < num.observed, num.observed, lower.ci)
upper.ci <- ans$n.hat + crit.value * ans$se.n.hat

ans$n.hat[1] <- NA  # can't estimate the first N
ans$se.n.hat[1] <- NA
lower.ci[1] <- NA
upper.ci[1] <- NA

# ----- Compute execution time here so can store in output.
ex.time <- (proc.time()[3] - start.sec) / 60

# ----- Done. Put into a 'CR' object.

ans <- list( histories=histories, 
    aux=aux, 
    intervals=intervals[1:(ns-1)], 
    loglike=ans$loglik, 
    deviance=ans$deviance, 
    aic=aic, 
    qaic=qaic, 
    aicc=aicc, 
    qaicc=qaicc,
    vif=ans$vif, 
    chisq.vif=ans$chisq.vif, 
    vif.df=ans$df.vif, 
    parameters=ans$parameters, 
    se.param=ans$se.param, 
    capcoef=capcoef, 
    se.capcoef=se.capcoef, 
    surcoef=surcoef, 
    se.surcoef=se.surcoef, 
    covariance=covariance,
    p.hat=ans$p.hat, 
    se.p.hat=ans$se.p.hat, 
    s.hat=ans$s.hat, 
    se.s.hat=ans$se.s.hat, 
    df=df, 
    df.estimated=ans$df.estimated,
    control=control,
    message=c(alg.mess,exit.mess,cov.mess), 
    exit.code=ans$exit.code, 
    cov.code=ans$cov.code, 
    fn.evals=ans$maxfn,
    ex.time = ex.time,
    n.hat =ans$n.hat, 
    se.n.hat=ans$se.n.hat, 
    n.hat.lower = lower.ci, 
    n.hat.upper = ceiling(upper.ci),
    n.hat.conf = conf, 
    nhat.v.meth = ans$nhat.v.meth, 
    num.caught=num.observed)
    
class(ans) <- c("cjs", "cr")



#   Compute fitted and residual components, if converged
if( ans$exit.code == 1 ){
    ans$fitted <- predict( ans )
    ans$residuals <- residuals( ans, type="pearson" )  
    ans$resid.type <- "pearson"
} else {
    ans$fitted <- NA
    ans$residuals <- NA  
    ans$resid.type <- NA
}

if( control$trace ){
    if( ex.time < 0.01666667 ){
        cat("\t(Execution time < 1 second)\n")
    } else {
        cat(paste("\t(Execution time =", round(ex.time,2), "minutes)\n"))
    }
}


ans

}
