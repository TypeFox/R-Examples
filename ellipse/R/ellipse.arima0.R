ellipse.arima0<-
  function(x, which = c(1, 2), level = 0.95, t = sqrt(qchisq(level, 2)), ...)
{
        ellipse.default(x$var.coef[which, which], centre = x$coef[which], t = t, ...)
}

