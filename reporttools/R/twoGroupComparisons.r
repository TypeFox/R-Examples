twoGroupComparisons <- function(vars, v0, conf.level = 0.95, paired = FALSE){

vars <- correctVarNames(vars, cols = NA)
tab <- data.frame(matrix(NA, nrow = ncol(vars), ncol = 16))
nam1 <- c("$n$", "$\\bar{x}$", "$\\tilde{x}$")
v0 <- as.factor(v0)
colnames(tab) <- c("Variable", paste(nam1, "$_{\\text{", levels(v0)[1], "}}$", sep = ""), paste(nam1, "$_{\\text{", levels(v0)[2], "}}$", sep = ""), "$n$ complete both", "mean difference", "Wilcoxon effect", "$p$-value $t$-test", 
    "$p$-value Wilcoxon", "CI $t$ low", "CI $t$ up", "CI Wilcoxon low", "CI Wilcoxon up")


for (i in 1:ncol(vars)){
    v <- vars[, i]
    
    v1 <- v[v0 == levels(v0)[1]]
    v2 <- v[v0 == levels(v0)[2]]    

    ## t-test
    t1 <- t.test(v1, v2, conf.int = TRUE, conf.level = conf.level, paired = paired)

    ## wilcoxon
    t2 <- try(wilcox.test(v1, v2, conf.int = TRUE, conf.level = conf.level, paired = paired), silent = TRUE)
    
    tab[i, 1]  <- colnames(vars)[i]
    tab[i, 2]  <- sum(is.na(v1) == FALSE)
    tab[i, 3]  <- mean(v1, na.rm = TRUE)
    tab[i, 4]  <- median(v1, na.rm = TRUE)
    tab[i, 5]  <- sum(is.na(v2) == FALSE)
    tab[i, 6]  <- mean(v2, na.rm = TRUE) 
    tab[i, 7]  <- median(v2, na.rm = TRUE)
    if (identical(paired, FALSE)){tab[i, 9]  <- -diff(t1$estimate)}
    if (identical(paired, TRUE)){
        tab[i, 8]  <- sum(complete.cases(cbind(v1, v2)))
        tab[i, 9]  <- t1$estimate
        }
    tab[i, 10]  <- t2$estimate
    tab[i, 11] <- t1$p.value
    tab[i, 12] <- t2$p.value
    tab[i, 13:14] <- t1$conf.int
    tab[i, 15:16] <- t2$conf.int
}

# format table
tab2 <- data.frame(matrix(NA, nrow = ncol(vars), ncol = 14))
colnames(tab2)[1:12] <- colnames(tab)[1:12]
colnames(tab2)[13:14] <- c("CI $t$", "CI Wilcoxon")

tab2[, 1] <- tab[, 1]
tab2[, c(2, 5, 8)] <- apply(tab[, c(2, 5, 8)], 1:2, disp, 0)
tab2[, c(11, 12)] <- apply(tab[, c(11, 12)], 1:2, formatPval, 2)
tab2[, c(3:4, 6:7, 9:10)] <- apply(tab[, c(3:4, 6:7, 9:10)], 1:2, disp, 2)
tab2[, 13] <- displayCI(as.matrix(tab[, 13:14]), 2)
tab2[, 14] <- displayCI(as.matrix(tab[, 15:16]), 2)

res <- list("raw" = tab, "formatted" = tab2)
return(res)
}















#
