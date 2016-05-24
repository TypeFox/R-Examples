pamer.fnc<-function (model, ndigits = 4) 
{
  if (length(rownames(anova(model))) == 0) {
    cat("nothing to evaluate: model has only an intercept.\n\n")
    cat("printing model fixed effects:\n")
    fixef(model)
  }
  else {
    dims <- NULL
    rank.X = qr(model@pp$X)$rank
    anova.table = anova(model)
    anova.table = cbind(anova.table, upper.den.df = nrow(model@frame) - rank.X)
    p.values.upper = as.numeric()
    p.values.lower = as.numeric()
    for (term in row.names(anova.table)) {
      p.values.upper = c(p.values.upper, round(1 - pf(anova.table[term,"F value"], anova.table[term, "Df"], nrow(model@frame) - rank.X), ndigits))
      model.ranef <- ranef(model)
      lower.bound <- 0
      for (i in 1:length(names(model.ranef))) {
        dims <- dim(model.ranef[[i]])
        lower.bound <- lower.bound + dims[1] * dims[2]
      }
      p.values.lower = c(p.values.lower, 1 - pf(anova.table[term, "F value"], anova.table[term, "Df"], nrow(model@frame) - rank.X - lower.bound))
    }
    dv <- gsub(" ", "", gsub("(.*)~.*", "\\1", as.character(model@call)[2]))
    ss.tot <- sum((model@frame[, dv] - mean(model@frame[, dv]))^2)
    aov.table <- as.data.frame(anova(model))
    expl.dev <- vector("numeric")
    for (i in rownames(aov.table)) {
      expl.dev <- c(expl.dev, aov.table[i, 2]/ss.tot)
    }
    names(expl.dev) <- rownames(aov.table)
    anova.table = round(cbind(anova.table, upper.p.val = p.values.upper, lower.den.df = nrow(model@frame) - rank.X - lower.bound, lower.p.val = p.values.lower, `expl.dev.(%)` = expl.dev * 100), ndigits)
    return(anova.table)
  }
}
