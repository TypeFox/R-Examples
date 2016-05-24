lassocoef <- 
function (formula, data, sopt, plot.opt = TRUE, ...) 
{
#    require(pls)
#    require(lars)
    mf <<- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- delintercept(model.matrix(mt, mf))

# ---- comment BL 14.2.2012 ---    
    mod_lasso <- lars(X, y)                                                     ## calculates entire lasso path
    aa <- apply(abs(mod_lasso$beta), 1, sum)                                    ## sum of absolute beta_i for i=1,...,m
    ind <- which.min(abs(aa/max(aa) - sopt))                                    ## aa/max(aa) == calculate "fraction" value; then substract sopt (given in fractions)
                                                                                ### reason: find index of sopt -> trick: search for index where (aa/max(aa) - sopt) = 0
                                                                                 ### thus take absolute values and search for minimum      
     coef <- predict(mod_lasso, s=sopt, type="coefficients", mode="fraction")$coefficients
    
    numb.zero <- sum(coef == 0)
    numb.nonzero <- sum(coef != 0)
    if (plot.opt) {
        plot(mod_lasso, breaks = FALSE, cex = 0.4, col = gray(0.6), 
            ...)
        abline(v = sopt, lty = 2)
        title(paste(numb.zero, "coefficients are zero,", numb.nonzero, 
            "are not zero"))
    }
    list( coefficients = coef,
          sopt = sopt, 
          numb.zero = numb.zero, 
          numb.nonzero = numb.nonzero, 
          ind = ind)
}

