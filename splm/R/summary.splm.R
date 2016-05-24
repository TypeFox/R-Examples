
`summary.splm` <-
function(object,...){

  ## summary method for splm objects
  ## adds incrementally to the model object, as summary.plm does
  ## structure remains the same for all type but 'spsegm' (symultaneous equations requires a special printing)
            ## to date, only balanced panels are allowed for 'splm'
            balanced <- TRUE #attr(object,"pdim")$balanced
            model.name <- object$type #attr(object,"pmodel")$model
            effect <- "individual" #attr(object,"pmodel")$effect
            ## make coefficients' table if vcov exist
            if (!is.null(object$vcov)) {
                std.err <- sqrt(diag(object$vcov))
               
#if(object$type == "fixed effects sarar")  std.err <- c(object$se.spat, sqrt(diag(object$vcov)))
                 #vcov(object) doesn't work
                b <- coefficients(object)
                z <- b/std.err
                p <- 2*pnorm(abs(z),lower.tail=FALSE)
                CoefTable <- cbind(b,std.err,z,p)
                colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
                object$CoefTable <- CoefTable
            } 
            else {
                object$CoefTable <- cbind(coefficients(object))
                colnames(object$CoefTable) <- c("Estimate")
            }

            # if (object$type == "fixed effects error" && object$method != "eigen") {
                # lambda <- object$spat.coef
                # object$lambda <- lambda
            # }

            if (object$type == "random effects GM" ) {
                lambda <- object$rho
                object$lambda <- lambda
            }

            if (object$type == "fixed effects GM" ) {
                lambda <- object$rho
                object$lambda <- lambda
            }

            ## make AR coefficient of y's table
            if(!is.null(object$vcov.arcoef)) {
                std.err1 <- sqrt(diag(object$vcov.arcoef))
                ar <- object$arcoef
                z <- ar/std.err1
                p <- 2*pnorm(abs(z),lower.tail=FALSE)
                ARCoefTable <- cbind(ar,std.err1,z,p)
                colnames(ARCoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
                object$ARCoefTable <- ARCoefTable
            }


            ## make error comps' table
            if(!is.null(object$vcov.errcomp)) {
                std.err2 <- sqrt(diag(object$vcov.errcomp))
                ec <- object$errcomp
                z <- ec/std.err2
                p <- 2*pnorm(abs(z),lower.tail=FALSE)
                ErrCompTable <- cbind(ec,std.err2,z,p)
                colnames(ErrCompTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
                object$ErrCompTable <- ErrCompTable
            }

            object$ssr <- sum(residuals(object)^2)
            object$tss <- tss(object$model[[1]])
            object$rsqr <- 1-object$ssr/object$tss
            object$fstatistic <- "nil" #Ftest(object)
            class(object) <- c("summary.splm","splm")
            object
        
}

