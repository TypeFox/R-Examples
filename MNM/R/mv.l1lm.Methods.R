`print.mvl1lm` <-
function(x, digits = 3, ...)
    {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)

    }


    
`coef.mvl1lm` <-
function(object,...)
    {
    object$coefficients
    }

`fitted.mvl1lm` <-
function(object,...)
    {
    xx <- object$fitted.values
    napredict(object$na.action, xx)
    }

`residuals.mvl1lm` <-
function(object,...)
    {
    xx <- object$residuals
    naresid(object$na.action, xx)
    }

`vcov.mvl1lm` <-
function(object,...)
    {
    object$vcov
    }

`predict.mvl1lm` <- function(object, newdata, na.action = na.pass, ...)
    {
    if (missing(newdata)) 
        return(object$fitted.values)
    
    
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
         xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    
    if (object$scores !="rank"){
         pred <- X %*% object$coefficients
          } else {
          if (object$IntC ==FALSE){
               pred <- X %*% object$coefficients
               } else {
                    coefs <- rbind(object$intercept, object$coefficients)
                    pred <- X %*% coefs
                    }
               }
               
    pred
    }


`summary.mvl1lm` <-
function(object,..., digits=3)
    {  
    fp <- format.pval(object$p.value, digits = 4)
    COEF <- coef(object)
    ny <- ncol(COEF)
    ynames <- colnames(COEF)
    if (is.null(ynames)) {
        lhs <- object$terms[[2L]]
        if (mode(lhs) == "call" && lhs[[1L]] == "cbind") 
            ynames <- as.character(lhs)[-1L]
        else ynames <- paste("Y", seq_len(ny), sep = "")
    }
    ind <- ynames == ""
    if (any(ind)) 
        ynames[ind] <- paste("Y", seq_len(ny), sep = "")[ind]
    stderror <- matrix(sqrt(diag(object$vcov)), ncol = ny)        
    value <- vector("list", ny)
    names(value) <- paste("Response", ynames)
    for (i in seq(ny)) {
    x.mat <- cbind(COEF[, i],stderror[,i])
        rownames(x.mat) <- rownames(COEF)
        colnames(x.mat) <- c("Estimate", "Std. Error")
        value[[i]] <- x.mat
        }
    class(value) <- "listof"
    cat("\n")
    cat(object$method)
    cat("\n")
    cat("\nCall:\n", deparse(object$call), "\n", sep = "") 
    if(object$scores == "rank"){
         if (object$IntC){
              if (object$stand == "inner"){
                  cat("\nInner HL-estimator for the residuals (intercept):\n")
                  print(object$intercept, digits=digits)
                  # cat("\n")
                  }
              if (object$stand == "outer"){
                  cat("\nOuter HL-estimator for the residuals (intercept):\n")
                  print(object$intercept, digits=digits)}
                  # cat("\n")
                  }
              }
    if (length(coef(object))) 
        {
        cat("\nTesting that all coefficients = 0:\n")
        cat(paste("Q.2 = ", format(round(object$statistic,4))," with ", object$parameter, " df, p.value ", if (substr(fp, 1L, 1L) == 
        "<") fp else paste("=", fp), sep="")) 
        cat("\n\n")
        cat("Results by response:\n")
        cat("\n")
        print(value, digits=digits, right=TRUE)
        cat("\n")
        }
    else cat("No coefficients\n")
    cat("\n")
   
    invisible(object)
    }



`plot.mvl1lm` <-
function(x, captation="Residuals vs fitted",...)
    {
    X <- x$fitted.values
    Y <- x$residuals
    colnames(Y) <- paste("Residuals",colnames(Y))
    colnames(X) <- paste("Fitted",colnames(X))
    pairs2(X, Y, main=captation, ...)
    }
