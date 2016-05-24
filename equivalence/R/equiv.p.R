# $Id: equiv.p.R,v 1.3 2005/09/28 01:37:53 andrewr Exp $

"equiv.p" <-
  function(x, y, alpha=0.05) {
    lm.equiv <- lm(y ~ I(x - mean(x, na.rm=TRUE)))
    coef.tab <- coef(summary(lm.equiv))
    target <- c(mean(x, na.rm=TRUE), 1)
    t.quant <- qt(p=1-alpha/2, df=df.residual(lm.equiv))
    output <- as.list((coef.tab[,1] -
                       t.quant*(1 - 2 * (coef.tab[,1] > target)) * 
                       coef.tab[,2]) - target)
    names(output) <- c("Intercept", "Slope")
    return(output)
  }

