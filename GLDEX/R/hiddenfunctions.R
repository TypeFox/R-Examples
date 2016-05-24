.gl.parameter.tidy<-
function (lambda1, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl") 
{
if(length(lambda1)<2) {       
lambda1 <- c(lambda1, lambda2, lambda3, lambda4) }
    lambda1
}


.qdgl.fmkl<-
function (p, lambdas) 
{
     if (!gl.check.lambda.alt1(lambdas, param = "fmkl", vect = TRUE)) {
       return(rep(NA,length(p)))
    }

    outside.range <- !as.logical((p <= 1) * (p >= 0))
    u <- p[!outside.range]
    dens <- lambdas[2]/(p^(lambdas[3] - 1) + (1 - p)^(lambdas[4] - 
        1))
    dens
}



.qdgl.rs<-
function (p, lambdas) 
{
     if (!gl.check.lambda.alt1(lambdas, param = "rs", vect = TRUE)) {
       return(rep(NA,length(p)))
    }

    outside.range <- !as.logical((p <= 1) * (p >= 0))
    u <- p[!outside.range]
    dens <- lambdas[2]/(lambdas[3] * (p^(lambdas[3] - 1)) + lambdas[4] * 
        ((1 - p)^(lambdas[4] - 1)))
    dens
}

.qgl.rs<-
function (p, lambdas) 
{
    u <- p
    outside.range <- !as.logical((p <= 1) * (p >= 0))
    u <- p[!outside.range]
    lambda4 = lambdas[4]
    lambda3 = lambdas[3]
    lambda2 = lambdas[2]
    lambda1 = lambdas[1]
    quants <- lambda1 + (u^lambda3 - (1 - u)^lambda4)/lambda2
    result<-rep(NA,length(p))
    result[!outside.range]<-quants
    result
}

.qgl.fmkl<-
function (p, lambdas) 
{
    lambda4 = lambdas[4]
    lambda3 = lambdas[3]
    lambda2 = lambdas[2]
    lambda1 = lambdas[1]
    outside.range <- !as.logical((p <= 1) * (p >= 0))
    u <- p[!outside.range]
    if (lambda3 == 0) {
        if (lambda4 == 0) {
            quants <- lambda1 + (log(u) - log(1 - u))/lambda2
        }
        else {
            quants <- lambda1 + (log(u) - ((1 - u)^lambda4 - 
                1)/lambda4)/lambda2
        }
    }
    else {
        if (lambda4 == 0) {
            quants <- lambda1 + ((u^lambda3 - 1)/lambda3 - log(1 - 
                u))/lambda2
        }
        else {
            quants <- lambda1 + ((u^lambda3 - 1)/lambda3 - ((1 - 
                u)^lambda4 - 1)/lambda4)/lambda2
        }
    }
    result <- rep(NA,length(p))
    result[!outside.range] <- quants
    result
}

       
.First<-function() { .setGLDEXEnv(.runif.halton.seed = list())
    .setGLDEXEnv(.rnorm.halton.seed = list())
    .setGLDEXEnv(.runif.sobol.seed = list())
    .setGLDEXEnv(.rnorm.sobol.seed = list())}

.GLDEXEnv <- new.env(hash = TRUE)

.setGLDEXEnv <-
    function(...)
{
    x <- list(...)
    nm <- names(x)
     if (is.null(nm) || "" %in% nm)
        stop("all arguments must be named")
    sapply(nm, function(nm) assign(nm, x[[nm]],
                                 envir = .GLDEXEnv))
    invisible()
}

.getGLDEXEnv <-
    function(x = NULL, unset = "")
{
    if (is.null(x))
        x <- ls(all.names = TRUE, envir = .GLDEXEnv)
###     unlist(mget(x, envir = .GLDEXEnv, mode = "any",
###                 ifnotfound = as.list(unset)), recursive = FALSE)
    get(x, envir = .GLDEXEnv, mode = "any")
}



  