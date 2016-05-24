# $Id: residplot.R 1012 2006-11-14 22:25:06Z ggorjan $

residplot  <-  function(model, formula, ...)
  {
    data  <- expand.model.frame( model, formula, na.expand=TRUE)

    newform  <- eval(parse( text=paste("as.call(", "resid(model) ~",
                        formula[-1],")" )))

    plot( newform, data=data, ylab="Residuals")
    lines(lowess( newform, data=data ), col="red")
    bandplot(newform,data=data)
  }

