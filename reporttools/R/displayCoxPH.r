displayCoxPH <- function(mod, cap = "", lab = "mod", dig.coef = 2, dig.p = 1){

    cap0 <- paste("Cox regression model, $n = $ ", summary(mod)$n, ", number of events = ", summary(mod)$nevent, ".", sep = "")
    if (cap != ""){cap <- paste(cap0, " ", cap, sep = "")} else {cap <- cap0}
    
    tab1 <- data.frame(cbind(summary(mod)$conf.int), summary(mod)$coefficients)
    tab1 <- data.frame(cbind(apply(tab1[, c(5, 1)], 1:2, disp, dig.coef), displayCI(as.matrix(tab1[, c(3, 4)]), digit = dig.coef), formatPval(tab1[, 9], dig.p)))
    colnames(tab1) <- c("coef", "HR = exp(coef)", "95\\% CI", "$p$-value")

    tab1 <- correctVarNames(tab1, cols = NA)

    xtab2 <- xtable(tab1, caption = cap, label = lab, align = "lrrrr")
    print(xtab2, include.rownames = TRUE, include.colnames = TRUE, floating = TRUE, type = "latex", size = "footnotesize", table.placement = "h!", 
        sanitize.text.function = function(x){x})
    return(tab1)
}
