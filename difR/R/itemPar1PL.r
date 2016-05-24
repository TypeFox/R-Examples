library(ltm)
library(lme4)

itemPar1PL<-function (data, engine = "ltm", discr = 1) 
{
    if (engine != "ltm" & engine != "lme4") 
        return("Error: 'engine' must be either 'ltm' or 'lme4'!")
    else {
        J <- ncol(data)
        if (engine == "ltm") {
          if (!is.null(discr)) {
            const <- rbind(c(J + 1, 1))
            mod <- rasch(data, constraint = const)
}
else mod <- rasch(data)
            par <- summary(mod)$coefficients[1:J, 1:2]
        }
        else {
            N <- nrow(data)
            C <- ncol(data)
            y <- NULL
            for (i in 1:ncol(data)) y <- c(y, as.numeric(data[, 
                i]))
            person <- 1:N
            pp <- rep(person, C)
            pp <- as.factor(pp)
            items <- rep(1, N)
            for (i in 2:C) items <- c(items, rep(i, N))
            items <- as.factor(items)
            mod <- glmer(y ~ items + (1 | pp) - 1, family = binomial, 
                REML = FALSE)
            par <- summary(mod)@coefs[, 1:2]
            par[, 1] <- -par[, 1]
        }
        colnames(par) <- c("b", "se(b)")
        if (is.null(colnames(data)) == FALSE) 
            row <- colnames(data)
        else {
            row <- NULL
            for (i in 1:J) row <- c(row, paste("Item", i, sep = ""))
        }
        rownames(par) <- row
        return(par)
    }
}