bas.haz <- function(

        formula,
        haz.formula,
        data,
        coef,
        var,
        alpha=.05,
        FUN=exp
)
{
     p <- length(attr(terms(formula),"term.labels"))
   
     Call <- match.call()
     fit <- coxph(formula,data,iter=0,init=coef[1:p])
     fit$call$formula <- Call[[2]]
     fit$call$data <- Call[[4]]
     fit$call$init <- Call[[5]]

     ipd.survfit <- survfit(fit)
 
     meta.bas.haz <- mlma.study.bashaz(
                  haz.formula,
                  coef,
                  var,
                  ipd.survfit,
                  alpha=alpha,
                  FUN=FUN
)

return(list(ipd.survfit=ipd.survfit,meta.bas.haz=meta.bas.haz))
}


mlma.survfit <- function(
             formula,
             data,
             coef
)

{

     Call <- match.call()
     fit <- coxph(formula,data,iter=0,init=coef)
     fit$call$formula <- Call[[2]]
     fit$call$data <- Call[[3]]
     fit$call$init <- Call[[4]]

return(survfit(fit))
}

mlma.study.bashaz <- function(
                  haz.formula,
                  coef,
                  var,
                  survfit,
                  alpha=.05,
                  FUN=exp
)

{
    p <- dim(var)[1]
    X <- data.frame(time = survfit$time)
    names(X) <- all.vars(terms(haz.formula))
    X <- model.matrix(terms(haz.formula),X)
    
    C <- matrix(0,nrow(X),p)
    q <- ncol(X)

    C[,(p-q+1):p] <- X
    estimates <- ci(C,coef,var,alpha=alpha,f=FUN)

    estimates <- data.frame(
              est = estimates[2,],
              lower = estimates[1,],
              upper = estimates[3,],
              time = survfit$time
    )

return(estimates)
}

