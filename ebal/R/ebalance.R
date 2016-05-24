ebalance <-
function(
              Treatment,
              X,
              base.weight = NULL,
              norm.constant  = NULL,
              coefs = NULL ,
              max.iterations = 200,
              constraint.tolerance = 1,
              print.level=0
              ){

# Checks 
    if (sum(Treatment  != 1 & Treatment  != 0) > 0) {
        stop("Treatment indicator ('Treatment') must be a logical variable, TRUE (1) or FALSE (0)")
    }
    if (var(Treatment) == 0) {
        stop("Treatment indicator ('Treatment') must contain both treatment and control observations")
    }

    Treatment <- as.numeric(Treatment)
    X         <- as.matrix(X)
    
    if (sum(is.na(X))>0){
       stop("X contains missing data")
    }
    if (sum(is.na(Treatment))>0){
       stop("Treatment contains missing data")
    }

    if (length(Treatment) != nrow(X)) {
        stop("length(Treatment) != nrow(X)")
    }

    if (length(max.iterations) != 1 ) {
        stop("length(max.iterations) != 1")
    }
    if (length(constraint.tolerance) != 1 ) {
        stop("length(constraint.tolerance) != 1")
    }
    
# set up elements
ntreated  <- sum(Treatment==1)
ncontrols <- sum(Treatment==0)

    if (is.null(base.weight)) {
        base.weight = rep(1, ncontrols)
    }
    if ( length(base.weight) !=  ncontrols) {
        stop("length(base.weight) !=  number of controls  sum(Treatment==0)")
    }

co.x <- X[Treatment==0,]
co.x <- cbind(rep(1,ncontrols),co.x)

    if(qr(co.x)$rank != ncol(co.x)){
      stop("collinearity in covariate matrix for controls (remove collinear covariates)")
    }


tr.total <- apply(as.matrix(X[Treatment==1,]),2,sum)

    if (is.null(norm.constant)) {
        norm.constant <- ntreated
    }
    if (length(norm.constant) != 1) {
        stop("length(norm.constant) != 1")
    }

tr.total <- c(norm.constant,tr.total)

   if(is.null(coefs)) {
    coefs = c(log(tr.total[1]/sum(base.weight)),rep(0,(ncol(co.x)-1)))
   }
   
   if(length(coefs) != ncol(co.x)) {
    stop("coefs needs to have same length as number of covariates plus one")
   }
   
## run algo
eb.out <- eb(tr.total=tr.total,
               co.x=co.x,
               coefs=coefs,
               base.weight=base.weight,
               max.iterations=max.iterations,
               constraint.tolerance=constraint.tolerance,
               print.level=print.level
               )

 if(eb.out$converged == TRUE & print.level>=0) {
        cat("Converged within tolerance \n")
     }

z <- list(
          target.margins = tr.total,
          co.xdata = co.x,
          w=eb.out$Weights.ebal,
          coefs=eb.out$coefs,
          maxdiff=eb.out$maxdiff,
          norm.constant = norm.constant,
          constraint.tolerance=constraint.tolerance,
          max.iterations=max.iterations,
          base.weight=base.weight,
          print.level=print.level,
          converged=eb.out$converged
    )
        
class(z) <- "ebalance"
return(z)

}

