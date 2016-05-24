"monitor" <-
function (a, n.chains=dim(a)[2], trans=NULL, keep.all=FALSE, Rupper.keep=FALSE) {

## If keep.all=T:  a is a n x m x k array:
##   m sequences of length n, k variables measured
## If keep.all=F:  a is a 2n x m x k array (first half will be discarded)
##
## trans is a vector of length k:  "" if no transformation, or "log" or "logit"
## (If trans is not defined, it will be set to "log" for parameters that
## are all-positive and 0 otherwise.)
##
## If Rupper.keep=TRUE:  keep Rupper.  (Otherwise don't display it.)
  invlogit <- function (x) {1 / (1 + exp(-x))}
  nparams <- if(length(dim(a)) < 3) 1 else dim(a)[length(dim(a))]
  # Calculation and initialization of the required matrix "output"
  output <- matrix( , ncol = if(n.chains > 1){if(Rupper.keep) 10 else 9} else 7, nrow = nparams)
  if (length(dim(a))==2) a <- array (a, c(dim(a),1))
  if (!keep.all){
    n <- floor(dim(a)[1]/2)
    a <- a[(n+1):(2*n), , , drop = FALSE]
  }
  if (is.null(trans))
    trans <- ifelse ((apply (a<=0, 3, sum))==0, "log", "")
  for (i in 1:nparams){
    # Rupper.keep:  discard Rupper (nobody ever uses it)
    ai <- a[ , , i, drop = FALSE]
    if (trans[i]=="log"){
        conv.p <- conv.par(log(ai), n.chains, Rupper.keep=Rupper.keep) # reason????
        conv.p <- list(quantiles = exp(conv.p$quantiles),
            confshrink = conv.p$confshrink, n.eff = conv.p$n.eff)
    }
    else if (trans[i]=="logit"){
        if (!is.R()){
            logit <- function (x) { log(x /(1- x)) }
        }    
        conv.p <- conv.par(logit(ai), n.chains, Rupper.keep=Rupper.keep)
        conv.p <- list(quantiles = invlogit(conv.p$quantiles),
            confshrink = conv.p$confshrink, n.eff = conv.p$n.eff)
    }
    else conv.p <- conv.par(ai, n.chains, Rupper.keep=Rupper.keep)
    output[i, ] <- c(mean(ai), sd(as.vector(ai)),
      conv.p$quantiles, 
        if(n.chains > 1) conv.p$confshrink, 
        if(n.chains > 1) round(conv.p$n.eff, min(0, 1 - floor(log10(conv.p$n.eff))))
      )
  }
  if(n.chains > 1)
    dimnames(output) <- list(dimnames(a)[[3]], c("mean","sd",
        "2.5%","25%","50%","75%","97.5%", "Rhat", if(Rupper.keep) "Rupper","n.eff"))
  else
    dimnames(output) <- list(dimnames(a)[[3]], c("mean","sd",
        "2.5%","25%","50%","75%","97.5%"))
   return (output)
}
