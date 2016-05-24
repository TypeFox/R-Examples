GLD.lm.simu <-
function (formula, data, init.coeff, init.resid, param, maxit = 20000, 
    fun, method = "Nelder-Mead",optim.fun,fit=NULL, init=NULL) 
{   fit<<-fit
    init.mod <- lm(formula, data = data)
    rm(fit)
    value <- c(init.coeff, init.resid[2:4])
    y <- init.mod$model[, 1]
    x <- model.matrix(init.mod)
    check <- optim.fun(value, x, y, param)
    if ((is.na(check) | is.inf(check))) {
        repeat {
            init.mod <- lm(formula, data = data)
            y <- init.mod$model[, 1]
            x <- model.matrix(init.mod)
if(is.null(init)){
    init.1 <- init.mod$coeff
    init.2 <- fun(init.mod$resid)}

else if(!is.null(init) & length(init.mod$coeff)==length(init)){
    init.1 <- init
    init.2 <- fun(y-x%*%init)}

else if(!is.null(init) & length(init.mod$coeff)!=length(init)){
    init.len<-length(init)
    init.1 <- init[1:(init.len-4)]
    init.2 <- init[(init.len-3):init.len]}
    
            value <- c(init.1, init.2[2:4])
            check <- optim.fun(value, x, y, param)
            if (is.na(check) == FALSE & is.inf(check) == FALSE) {
                break
            }
        }
    }
        
    result <- optim(value, optim.fun, x = x, y = y, 
        param = param, control = list(maxit = 20000), method = method)
    full.result <- result$par[-c((length(result$par) - 2):length(result$par))]
    adj <- mean(y - x %*% full.result)
    if (is.element("(Intercept)", dimnames(x)[[2]]) == TRUE) {
        full.result[1] <- full.result[1] + adj
    }
    if (is.element("(Intercept)", dimnames(x)[[2]]) == FALSE) {
        warning("Bias adjustment is provided separately from estimated 
        parameters, adjustment is necessary to ensure residuals sum to zero")
    }
    full.result <- c(full.result, fun.mean.convert(c(0, 
    result$par[c((length(result$par) - 
        2):length(result$par))]), param))
    names(full.result)[1:(length(full.result) - 4)] <- dimnames(x)[[2]]
    names(full.result)[(length(full.result) - 3):length(full.result)] <- paste("L", 
        1:4, sep = "")
    if (result$convergence == 0) {
        converge.report <- "converged"
    }
    if (result$convergence != 0) {
        converge.report <- "not converged"
    }
    message1 <- paste(toupper(param), "GLD")
    message2 <- converge.report
    r <- cbind(adj, t(full.result), message1, message2)
    return(r)
}
