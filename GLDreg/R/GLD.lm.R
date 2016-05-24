GLD.lm<-
function (formula, data, param, maxit = 20000, fun, method = "Nelder-Mead", 
    diagnostics = TRUE, range = c(0.01, 0.99),init=NULL) 
{
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
    if (identical(fun, fun.RMFMKL.ml.m) | identical(fun, fun.RMFMKL.ml) | 
        identical(fun, fun.RPRS.ml.m) | identical(fun, fun.RPRS.ml)) {
        result <- optim(value, fun.model.lm.optim, x = x, y = y, 
            param = param, control = list(maxit = 20000), method = method)
        result <- optim(result$par, fun.model.lm.optim, x = x, 
            y = y, param = param, control = list(maxit = 20000), 
            method = method)
        msgr<-"Maximum Likelihood Estimation"
    }
    if (identical(fun, fun.RMFMKL.lm) | identical(fun, fun.RPRS.lm)) {
        result <- optim(value, fun.model.lm.optim.Lmoment, x = x, 
            y = y, param = param, control = list(maxit = 20000), 
            method = method)
        result <- optim(result$par, fun.model.lm.optim.Lmoment, 
            x = x, y = y, param = param, control = list(maxit = 20000), 
            method = method)
        msgr<-"L moment matching"
    }
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
    message1 <- paste("This analysis was carried out using", 
        toupper(param), "GLD")
    message2 <- paste("The error distribution was estimated using", 
                  msgr)
    message3 <- paste("The optimisation procedure used was", 
        as.character(substitute(method)), "and it has", converge.report)
    messages <- rbind(message1, message2, message3)
    dimnames(messages)[[1]] <- NULL
    fitted <- x %*% full.result[1:(length(full.result) - 4)]
    if (is.element("(Intercept)", dimnames(x)[[2]]) == TRUE) {
        resid <- y - fitted
    }
    if (is.element("(Intercept)", dimnames(x)[[2]]) == FALSE) {
        resid <- y - (fitted + adj)
    }
    if (diagnostics) {
        gld.values <- full.result[(length(full.result) - 3):length(full.result)]
        par(mfrow = c(1, 1))
        qqgld.default(resid, gld.values, param = param)
        legend("bottomright", c(paste(toupper(param), "GLD"), 
            paste("(", paste(signif(gld.values, 3), collapse = ","), 
                ")", sep = "")), bty = "n")
        pval <- ks.gof(resid, "pgl", lambda1 = gld.values, param = param)$p.value
        legend("topleft", paste("KS test\np-value=", format.pval(pval)), 
            bty = "n")
    }

print(messages)
print(full.result)

    return(list(Message = messages, `Bias Correction` = adj, 
        `Estimated parameters` = full.result, Fitted = fitted, 
        Residual = resid, formula = formula, param = param, y = y, 
        x = x, fun = fun))
}
