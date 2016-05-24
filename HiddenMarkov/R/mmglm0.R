"mmglm0" <-
function (x, Pi, delta, family, link, beta, glmformula=formula(y~x1),
          sigma=NA, nonstat=TRUE, msg=TRUE)
{
    if (msg)
        message("NOTE: 'mmglm1' is an updated and more general version of 'mmglm0'.")
    x <- c(list(x=x, Pi=Pi, delta=delta, family=family, link=link,
                beta=beta, glmformula=glmformula, sigma=sigma,
                nonstat=nonstat))
    class(x) <- c("mmglm0")
    return(x)
}


