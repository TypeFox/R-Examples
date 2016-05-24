comp.listEndorse <- function(y.endorse, y.list, treat, n.draws = 10000, alpha = .05, endorse.mean = FALSE,
                             method = "pearson") {

  if (endorse.mean == FALSE) {
    y.endorse.mean <- apply(y.endorse, 1, function(x) mean(x, na.omit = T))
  } else {
    y.endorse.mean <- y.endorse
  }
  
  calc.corr <- function(y.endorse.mean, y.list, treat) {
    
    no.na.vec <- !(is.na(y.endorse.mean) | is.na(y.list))
    
    cor.treat <- cor(y.endorse.mean[treat == 1 & no.na.vec == 1],
                     y.list[treat==1 & no.na.vec == 1], method = method)
    cor.control <- cor(y.endorse.mean[treat == 0 & no.na.vec == 1],
                       y.list[treat==0 & no.na.vec == 1], method = method)
    return(list(cor.treat = cor.treat, cor.control = cor.control))
  }

  calc.d <- function(cor.treat, cor.control) {
    
    d <- cor.treat - cor.control
    
    return(d)
  }

  cor <- calc.corr(y.endorse.mean, y.list, treat)

  d.bootstrap <- rep(NA, n.draws)
  for(i in 1:n.draws) {
    sample <- sample(1:length(treat), length(treat), replace = T)
    cor.bootstrap <- calc.corr(y.endorse.mean[sample], y.list[sample], treat[sample])
    d.bootstrap[i] <- calc.d(cor.treat = cor.bootstrap$cor.treat,
                               cor.control = cor.bootstrap$cor.control)
  }

  return(list(cor.treat = cor$cor.treat, cor.control = cor$cor.control,
              p.value = mean(d.bootstrap <= 0),
              ci = c(quantile(d.bootstrap, alpha/2), quantile(d.bootstrap, (1-alpha)/2))))
  
}


