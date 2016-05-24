
sphtest <- function (x, ...)
{
    UseMethod("sphtest")
}

sphtest.formula <- function (x, data, index = NULL, listw, spatial.model = c("lag", "error", "sarar"), method = c("ML", "GM"), errors = c("KKP", "BSK"),...)
{
 ## performs a Hausman test of a FE model with spatial lag or error
 ## against "alternative" with same spatial specification

    switch(match.arg(spatial.model),
    lag = {
    	lag = TRUE
    	spatial.error = FALSE
    	},
    error = {
    	lag = FALSE
    	spatial.error = TRUE
    	},
    sarar = {
    	lag = TRUE
    	spatial.error = TRUE
    	},
    	
    )

errors <- match.arg(errors)

    x0 <- update(x, .~.-1)

    method <- switch(match.arg(method), 

    ML = {
    	
    femod <- spml(x, data = data, index = index, listw = listw, lag = lag, spatial.error = spatial.error, model = "within", errors = errors)

    remod <- spml(x, data = data, index = index, listw = listw, lag = lag, spatial.error = spatial.error, model = "random", errors = errors)
  	
    	},
    	
    GM = {
    	
  femod <- spgm(x, data = data, index = index, listw = listw, lag = lag, spatial.error = spatial.error, model = "within", moments = "fullweights")
  
  remod <- spgm(x, data = data, index = index, listw = listw, lag = lag, spatial.error = spatial.error, model = "random", moments = "fullweights")
    	
    	},	
    
    stop("\n Unknown method")
    )    
    
    
    sphtest(femod, remod, ...)
}

sphtest.splm <- function (x, x2, ...){


  ## check that the models have the same specification but different effects
  
if (!all.equal(x$legacy, x2$legacy)) stop("The model are different")
if(x$ef.sph == x2$ef.sph) stop("Effects should be different")

    ran <- match("random", c(x$ef.sph, x2$ef.sph))

if(ran == 1){

	xwith <- x2
	xbetw <- x

	}    

if(ran == 2){

	xwith <- x
	xbetw <- x2

	}    	
    
  ## test on coefficients (excluding SAR)      
  ## model order is irrelevant
  
    tc <- match(names(coef(xwith)), names(coef(xbetw)) )

    coef.wi <- coef(xwith)
    coef.re <- coef(xbetw)[tc]
    vcov.wi <- xwith$vcov
    vcov.re <- xbetw$vcov[tc,tc]
    
    dbeta <- coef.wi - coef.re
    df <- length(dbeta)
    dvcov <- vcov.re - vcov.wi
    stat <- abs(t(dbeta) %*% solve(dvcov) %*% dbeta)
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    names(stat) <- "chisq"
    parameter <- df
    names(parameter) <- "df"
    data.name <- paste(deparse(x$call$formula))
    alternative <- "one model is inconsistent"
    res <- list(statistic = stat, p.value = pval, parameter = parameter,
        method = "Hausman test for spatial models", data.name = data.name, alternative = alternative)
    class(res) <- "htest"
    return(res)
}
