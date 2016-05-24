F.3d.model.matrix <- function( formula, d1, d2 ){
#
#   Produce one 3d model matrix for CR modeling
#
#   formula = formula object without or without response.
#       response is ignored.
#   d1 = magnitude of dimension 1 = number of rows = nan
#   d2 = Magnitude of dimension 2 = number of cols = ns

    call <- match.call()
    contrasts <- NULL
    mf <- match.call(expand.dots = FALSE)
    mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
    mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
    mf$... <- NULL
    mf$drop.unused.levels <- TRUE
    # Must wipe out names of extra arguments to this function
    mf$d1 <- NULL
    mf$d2 <- NULL
    mf$survival <- NULL
    mf[[1]] <- as.name("model.frame")
    mf$formula <- formula
    
    form <- as.character(formula)[2]
    if( nchar(form)==1 & form == "1" ){
	   # this is a model with intercept only
	   survX <- matrix(1,d1,d2)
       surv.intercept <- TRUE
	   ex.xvars <- "(Intercept)"
    } else {
       # this model could be empty, could have x variables, maybe an intercept
       mf <- eval( mf )
       mt <- attr(mf, "terms")
       surv.intercept <- attr( mt, "intercept") == 1
       if( (length(attr( mt, "order")) == 0) & !surv.intercept ){
            stop("Empty models not allowed")
       }
       xvars <- as.character(attr(mt, "variables"))[-1]
       xlev <- if (length(xvars) > 0) {
        	   xlev <- lapply(mf[xvars], levels)
       }
       xcontr <- lapply( mf[xvars], function(x){attr(x,"contr")} )

   	   survX <- NA
       survX <- model.matrix(mt, mf, contrasts)
    
       assign.col <- attr( survX, "assign" )
       assign.col <- assign.col[ assign.col > 0 ]  # don't want intercept
       n.cols <- as.numeric(table(assign.col))
       nv <- length( n.cols )  # number variables besides intercept
       ny <- 0
       ex.xvars <- NULL
       if( nv > 0 ){  # this handles cases like ~1 + x3 - x3
            for( i in 1:nv){
                if( n.cols[i] %% d2 != 0 ){ stop(paste(xvars[i], "expanded into wrong number of columns")) }
                ncoef <- n.cols[i] / d2
                if( ncoef == 1 ){
                    # not a factor
                    ex.xvars <- c(ex.xvars, xvars[i])
                    ny <- ny + 1
                } else {
                    # is a factor.  Levels to drop are taken care of in ivar and tvar
                    ex.xvars <- c(ex.xvars, paste(xvars[i], ":", xlev[[i]], sep=""))
                    ny <- ny + ncoef
                }
            }
       } 

       # Expand intercept into a matrix
       if(surv.intercept){
            survX <- cbind( matrix(1,nrow=d1,ncol=d2), survX[,-1] )
            ny <- ny + 1
            ex.xvars <- c("(Intercept)",ex.xvars)
       }

    }
    attr(survX, "intercept") <- surv.intercept
    attr(survX, "variables") <- ex.xvars

    survX
}
