multi.snp.test <-
  function(y, x, x.adj = NULL,
           type = c("gaussian", "binomial")) {    
    switch (type,
            gaussian =  test <- "F",
            binomial =  test <- "Chisq")
    not.adjusted <- FALSE
    if (all(is.null(x.adj))) {
      not.adjusted <- TRUE
    }
    if (not.adjusted) {
      # no covariates (only intercept)
      fit.glm0 <- glm(y ~ 1, family = type)
      fit.glm1 <- glm(y ~ ., data = as.data.frame(x), family = type)
      aov.glm  <- anova(fit.glm0, fit.glm1, test = test)
    } else {
      # covariates
      if (!is.data.frame (x.adj))
        x.adj <- as.data.frame (x.adj, stringsAsFactors = FALSE)
      if (nrow(x) != nrow(x.adj)) {
        stop (paste("Don't match nrow(x) = ", 
                    nrow(x),
                    " and nrow(x.adj) = ",
                    nrow(x.adj), 
                    collapse = ""))
      }
      x.all <- cbind (x, x.adj)
      fit.glm0 <- glm (y ~ ., data = x.adj, family = type)
      fit.glm1 <- glm (y ~ ., data = x.all, family = type)
      aov.glm  <- anova (fit.glm0, fit.glm1, test = test)
    }
    return(list(fit.glm0 = fit.glm0,
                fit.glm1 = fit.glm1,
                aov.glm  = aov.glm,
                beta     = coef(summary(fit.glm1))))
  }
