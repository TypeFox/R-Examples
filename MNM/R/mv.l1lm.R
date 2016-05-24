#################################################################################################################################################
#
# Function for regression calling subfunctions
#
#################################################################################################################################################


mv.l1lm <- function(formula, scores = "identity", stand = "outer", maxiter = 1000, eps = 1e-6, eps.S = 1e-6, x = TRUE, y = TRUE, data, subset, na.action)
{
    ret.x <- x
    ret.y <- y

    oldC <- getOption("contrasts")
    options(contrasts = c("contr.treatment", "contr.treatment"))
    cl <- match.call()
    ## Get the model frame
    mf <- match.call(expand.dots = FALSE)
    
    scores <- match.arg(scores, c("identity", "sign", "rank"))
    stand <- match.arg(stand, c("outer", "inner"))
    
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    
    
    ## Get the terms
    mt <- attr(mf, "terms")        # This is to include the `predvars'
                                   # attribute of the terms
    ## Get the data matrices
    Y <- model.response(mf, "numeric")
    if (is.vector(Y)) stop("response must be at least bivariate")
    if (is.matrix(Y)) {
        if (is.null(colnames(Y)))
            colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    } else {
        Y <- as.matrix(Y)
        colnames(Y) <- deparse(formula[[2]])
    }
    X <- model.matrix(mt, mf)


res1<-switch(scores,
        "identity"={
               reg.identity(Y=Y, X=X)
               }
        ,
        "sign"={
               switch(stand,
                    "outer" = {
                    reg.signs.outer(Y=Y, X=X, eps=eps, maxiter=maxiter, eps.S=eps.S)
                    }
               ,
                    "inner" = {
                    reg.signs.inner(Y=Y, X=X, eps=eps, maxiter=maxiter, eps.S=eps.S)
                    }
                    )
                    }
        ,
        "rank"={
               switch(stand,
                    "outer" = {
                    reg.ranks.outer(Y=Y, X=X, eps=eps, maxiter=maxiter, eps.S=eps.S)
                    }
               ,
                    "inner" = {
                    reg.ranks.inner(Y=Y, X=X, eps=eps, maxiter=maxiter, eps.S=eps.S)
                    }
                    )
                    }
        )
        options(contrasts = oldC)

CheckR <- rowSums(res1$residuals^2) < eps
if (sum(CheckR) > ncol(Y)+1) warning(paste(sum(CheckR), " > (p+1) residuals have a smaller norm than ", eps,". This might be a problem!", sep="")) 

res<-c(res1, list(call=cl, na.action = attr(mf, "na.action"), xlevels = .getXlevels(mt, mf), terms = mt, contrasts <- attr(X, "contrasts")))
if (ret.x)  res <- c(res, list(x = X)) else res <- list(res, list(x = NULL))
if (ret.y)  res <- c(res, list(y = Y)) else res <- list(res, list(y = NULL))
class(res) <- "mvl1lm"
return(res)
}
