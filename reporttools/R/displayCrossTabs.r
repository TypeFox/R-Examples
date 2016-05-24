displayCrossTabs <- function(vars, v0, nam0, lab0, percentage = c("none", "row", "col", "total")[1], add.p = TRUE){

tabs <- NULL
ps <- NULL
tests <- NULL
vars <- correctVarNames(vars, cols = NA)

for (i in 1:ncol(vars)){

    vs <- data.frame(v0, vars[i])
    vs <- eliminateNA(vs)$complete
    v1 <- factor(vs[, 1], exclude = NULL)
    v2 <- factor(vs[, 2], exclude = NULL)
    capo <- paste(nam0, " vs. ", colnames(vars)[i], ".", sep = "")
    p1 <- NA
    fish <- NA
    
    if (identical(add.p, TRUE)){
        if ((nlevels(v1) > 1) & (nlevels(v2) > 1)){
            fish <- "$\\chi^2$"
            chi <- suppressWarnings(chisq.test(v1, v2))
            p1 <- chi$p.value
            if (min(chi$expected) <= 5){
                fish <- "Fisher's exact"
                p1 <- fisher.test(v1, v2)$p.value
                }
            capo <- paste(nam0, " vs. ", colnames(vars)[i], ".", sep = "")
            capo <- paste(capo, " $p$-value ", fish, " test: ", disp(p1), ".", sep = "")
            }
    }
    displayKbyC(v1, v2, names = c(nam0, colnames(vars)[i]), cap = capo, lab = paste(lab0, i, sep = ""), percentage = percentage)
    tabs[[i]] <- table(v1, v2)
    ps[i] <- p1
    tests[i] <- fish
    }
    
    res <- list("tabs" = tabs, "ps" = ps, "tests" = tests)
    return(res)
}
