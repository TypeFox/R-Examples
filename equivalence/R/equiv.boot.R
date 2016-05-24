# $Id: equiv.boot.R,v 1.3 2005/10/10 10:14:43 andrewr Exp $

"equiv.boot" <-
  function(x, y, alpha=0.05, b0.ii=0.25, b1.ii=0.25, reps=100,
            b0.ii.absolute = FALSE) {
    if (!(length(x) == length(y))) stop ("Data must be paired.") 
    plot.boot <- boot(cbind(y, x),
                      equiv.boot.lm,
                      R=reps,
                      rel.int.int=b0.ii, rel.int.slope=b1.ii,
                    b0.ii.absolute = b0.ii.absolute)
    eff.alpha <- 1 - sqrt(1-alpha)
    c.b0.l <- quantile(plot.boot$t[,7], probs=eff.alpha)
    c.b0.u <- quantile(plot.boot$t[,7], probs=1-eff.alpha)
    i.b0.l <- ifelse(b0.ii.absolute,
                     mean(x, na.rm = TRUE) - b0.ii,
                     mean(x, na.rm = TRUE) * (1 - b0.ii))
    i.b0.u <- ifelse(b0.ii.absolute,
                     mean(x, na.rm = TRUE) + b0.ii,
                     mean(x, na.rm = TRUE) * (1 + b0.ii))
    c.b1.l <- quantile(plot.boot$t[,8], probs=eff.alpha)
    c.b1.u <- quantile(plot.boot$t[,8], probs=1-eff.alpha)
    i.b1.l <- 1 * (1 - b1.ii)
    i.b1.u <- 1 * (1 + b1.ii)
    quantiles <- apply(plot.boot$t, 2, mean)[1:6]
    list(n = sum(complete.cases(cbind(x, y))),
         ci.b0 = c(c.b0.l, c.b0.u),
         rs.b0 = c(i.b0.l, i.b0.u),
         q.b0 = quantiles[1:3],
         Test.b0 = ifelse(c.b0.l > i.b0.l & c.b0.u < i.b0.u,
           "Reject", "not Reject"),
         ci.b1 = c(c.b1.l, c.b1.u),
         rs.b1 = c(i.b1.l, i.b1.u),
         q.b1 = quantiles[4:6],
         Test.b1 = ifelse(c.b1.l > i.b1.l & c.b1.u < i.b1.u,
           "Reject", "not Reject"),
         eff.alpha = eff.alpha
         )
  }


