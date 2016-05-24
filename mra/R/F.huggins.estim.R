F.huggins.estim <- function(capture, recapture=NULL, histories, remove=FALSE, cap.init, recap.init,
    nhat.v.meth=1, df=NA, link="logit", control=mra.control()){

start.sec <- proc.time()[3]

run.date <- Sys.time()

if( missing(histories) ){
    stop("Capture histories matrix must be specified")
}
if( missing(capture) ){
    stop("Capture covariates must be specified")
}

if( length(union( unique(histories), c(0,1))) > 2 ) stop("Capture histories must consist of 0's and 1's only.")

# Remove rows of all zeros. These are errors, and we could stop, but I'll just remove.
zero.ind <- apply( histories, 1, sum ) == 0
if( any(zero.ind) ){
    stop(paste("Rows of all zeros not allowed in history matrix.", sum(zero.ind), "rows found."))
}


hist.name <- deparse(substitute(histories))
cr.call <- match.call()
nan <- nrow( histories )
ns <- ncol( histories )

# return the X and Y matrix.  Model for initial captures is in 'capture'. 
# Model for recaptures is in 'recapture'. 
# After this, both matrices are huge NAN by (ns*nx) matrices.

capX <- F.3d.model.matrix( as.formula(capture), nan, ns )
cap.intercept <- attr(capX, "intercept") == 1
cap.names <- attr(capX, "variables")
nx <- length(cap.names)

if( !is.null(recapture) ){
    recapX <- F.3d.model.matrix( as.formula(recapture), nan, ns )
    recap.intercept <- attr(recapX, "intercept") == 1
    recap.names <- attr(recapX, "variables")
    ny <- length(recap.names)
} else {
    recapX <- rep(1,nan)
    recap.intercept <- FALSE
    recap.names <- "NULL"
    ny <- 0
}

# Rep out the remove vector, and convert to integer
if(length(remove) < nx){
    remove <- rep(remove, nx)
} else if (length(remove) > nx ){
    remove <- remove[1:nx]
}
remove.vec <- as.numeric( remove )

#   Set initial values if missing or short
if( missing(cap.init) ){
    cap.init <- rep(0,nx)
} else if(length(cap.init) < nx ){
    cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 

if( missing(recap.init) ){
    recap.init <- rep(0,ny)
} else if(length(recap.init) < ny ){
    recap.init <- c(recap.init, rep(0, ny-length(recap.init)))
} 

#   Set up the tolerance vector, if not specified, or if not long enough
if( length(control$tol) < (nx+ny) ){
    control$tol <- rep(control$tol, trunc((nx+ny) / length(control$tol))+1)[1:(nx+ny)]
} else if( length(control$tol > (nx+ny)) ){
    control$tol <- control$tol[1:(nx+ny)]
}



#   Do the estimation, but first allocate room for answers
loglik <- deviance <- aic <- qaic  <- lower.ci <- upper.ci <- 0
parameters <- se.param  <- rep(0, nx+ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- c.hat <- se.c.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- 0
if( is.na(df) ){
    df.estimated <- 1  # Have MRA estimate rank of var-covar matrix
} else {
    df.estimated <- 0  # Don't bother, df either set by user or will use nx+ny
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


if( control$trace ) cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")

ans <- .Fortran( "hugginsmodel", 
        nan         = as.integer(nan), 
        ns          = as.integer(ns), 
        nx          = as.integer(nx), 
        ny          = as.integer(ny),
        histories   = as.integer(histories),  
        remove.vec  = as.integer(remove.vec),
        algorithm   = as.integer(control$algorithm), 
        cov.meth    = as.integer(control$cov.meth), 
        nhat.v.meth = as.integer(nhat.v.meth), 
        capX        = as.double(capX),
        recapX      = as.double(recapX), 
        cap.init    = as.double(cap.init),
        recap.init  = as.double(recap.init), 
        trace       = as.integer(control$trace),
        link        = as.integer(link.code),
        maxfn       = as.integer(control$maxfn),
        tol         = as.double(control$tol),
        loglik      = as.double(loglik), 
        deviance    = as.double(deviance), 
        aic         = as.double(aic),   
        parameters  = as.double(parameters),
        se.param    = as.double(se.param), 
        covariance  = as.double(covariance), 
        p.hat       = as.double(p.hat), 
        se.p.hat    = as.double(se.p.hat), 
        c.hat       = as.double(c.hat), 
        se.c.hat    = as.double(se.c.hat),
        n.hat       = as.double(n.hat), 
        se.n.hat    = as.double(se.n.hat), 
        lower.ci    = as.double(lower.ci),
        upper.ci    = as.double(upper.ci),
        exit.code   = as.integer(exit.code), 
        cov.code    = as.integer(cov.code), 
        df.estimated= as.integer(df.estimated), 
        PACKAGE="mra" )

if(control$trace) cat(paste("Returned from MRA. Details in MRA.LOG.\n", sep=""))

# ----- Fortran sets missing standard errors < 0. Reset missing standard errors to NA.
ans$se.param[ ans$se.param < 0 ] <- NA

# ----- R does not preserve the matrix structures in .Fortran call.  Put matricies, 
#   which are now vectors, back to matricies.
ans$covariance <- matrix( ans$covariance, nrow=nx+ny ) 
ans$p.hat      <- matrix( ans$p.hat, nrow=nan )
ans$se.p.hat   <- matrix( ans$se.p.hat, nrow=nan )
ans$c.hat      <- matrix( ans$c.hat, nrow=nan )
ans$se.c.hat   <- matrix( ans$se.c.hat, nrow=nan )



# ----- Work out exit codes.  Code is returned by VA09AD.
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

if(ans$cov.meth == 1){
    cov.mess = "Covariance from numeric derivatives."
} else if (ans$cov.meth == 2){
    cov.mess = "Covariance from optimization Hessian."
} else {
    cov.mess = "Unkown covariance method."
}

if(ans$cov.code == 0){
    cov.mess = paste( cov.mess,  "SUCCESS: Non-singular covariance matrix.")
} else if( cov.code == 1) {
    cov.mess = paste( cov.mess,  "ERROR: COVARIANCE MATRIX IS SINGULAR.")
} else {
    cov.mess = paste( cov.mess,  "ERROR: NON-NEGATIVE DEFINITE COVARIANCE MATRIX.")
}

# ----- Remove trailing blanks from message
message <- c( alg.mess, exit.mess, cov.mess)
if( control$trace ) cat(paste("\t(", message, ")\n", sep=""))

# ----- Wipe out the first recapture probability.  Can't have recapture in first period.
ans$c.hat[,1] <- NA
ans$se.c.hat[,1] <- NA

# ----- Transfer over coefficients
capcoef <- ans$parameters[1:nx]
se.capcoef <- ans$se.param[1:nx]
names(capcoef) <- cap.names
names(se.capcoef) <- cap.names
nms <- paste("cap.",names( capcoef ),sep="")

if( ny >= 1 ){
    recapcoef <- ans$parameters[(nx+1):(nx+ny)]
    se.recapcoef <- ans$se.param[(nx+1):(nx+ny)]
    names(recapcoef) <- recap.names
    names(se.recapcoef) <- recap.names
    nms <- c( nms, paste("recap.",names( recapcoef ), sep="") )
} else {
    recapcoef <- NULL
    se.recapcoef <- NULL
}    

dimnames(covariance) <- list( nms, nms )

# ----- Fix up number of parameters.
if( is.na(df) ){
    df <- ans$df.estimated   # use the rank estimated by MRAWIN
} else if( df <= 0 ){
    df <- nx+ny  # assume full rank
} # else use the values supplied by user (unchanged from input)

# ----- Now that df is computed, recompute fit statistics
aic <- -2*ans$loglik + 2*df
n.eff <- nan * ns
aicc <- aic + (2*df*(df+1))/(n.eff - df - 1)



# ----- Put all the "auxillary" info together
aux <- list( call=cr.call, nan=nan, ns=ns, nx=length(capcoef), ny=length(recapcoef),
        cov.name=nms, 
        ic.name=hist.name, 
        mra.version=packageDescription("mra")$Version, 
        R.version=R.version.string,
        run.date=run.date )


# ----- Estimates of N are computed in Fortran.  No need to modify.

# ----- Done. Put into a 'hug' object.
ex.time <- (proc.time()[3] - start.sec) / 60
ans <- list( histories=histories, 
    aux=aux, 
    loglike=ans$loglik, 
    deviance=ans$deviance, 
    aic=aic, 
    aicc=aicc, 
    capcoef=capcoef, 
    se.capcoef=se.capcoef, 
    recapcoef=recapcoef, 
    se.recapcoef=se.recapcoef,
    remove=remove,
    covariance=ans$covariance,
    p.hat=ans$p.hat, 
    se.p.hat=ans$se.p.hat, 
    c.hat=ans$c.hat, 
    se.c.hat=ans$se.c.hat,
    df=df, 
    control=control,
    message=message, 
    fn.evals = ans$maxfn,
    ex.time = ex.time,
    exit.code=ans$exit.code, 
    cov.code=ans$cov.code, 
    cov.meth=ans$cov.meth, 
    n.hat = ans$n.hat, 
    se.n.hat = ans$se.n.hat, 
    n.hat.lower = ans$lower.ci, 
    n.hat.upper = ans$upper.ci,
    n.hat.conf = 0.95, 
    nhat.v.meth = ans$nhat.v.meth, 
    num.caught=nan,
    n.effective=n.eff)
class(ans) <- c("hug", "cr")

#need to check this returned object with the documentation.  

if( control$trace ) {
    if( ex.time < 0.01666667 ){
        cat("\t(Execution time < 1 second)\n")
    } else {
        cat(paste("\t(Execution time =", round(ex.time,2), "minutes)\n"))
    }
}

ans

}
