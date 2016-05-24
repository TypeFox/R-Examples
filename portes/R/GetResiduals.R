"GetResiduals" <-
  function(obj)
{
    if (class(obj) != "ar" && class(obj) != "arima0" && class(obj) != "Arima" && class(obj) != "varest" && class(obj) != "FitAR" 
       && class(obj) != "FitFGN" && class(obj) != "garch" && class(obj) != "fGARCH" && class(obj) != "list" ) 
        stop("must be class ar, arima0, Arima, varest, FitAR, FitFGN, garch, fGARCH, or list object")
    if (class(obj)=="ar"){
        order <- obj$order
        res <- ts(as.matrix(obj$resid)[-(1:order),])
    }
    else if (all(class(obj) == "arima0") || all(class(obj) == "Arima")) {
	  pdq <- obj$arma
	  p <- pdq[1]
	  q <- pdq[2]
	  d <- pdq[6]
        order <- p+q
         res <- ts(obj$residuals) 
    }
    else if (class(obj)=="varest"){
     order <- obj$p
     res <- resid(obj)
    }
    else if (class(obj)=="FitAR"){ 
     order <- length(obj$phiHat)
     res <- ts(obj$res) 
    }
    else if (class(obj) == "FitFGN") {
        order <- 0
        res <- ts(obj$res)
    }
    else if (class(obj) == "garch") {
	  GarchOrder <- as.vector(obj$order)
        order <- 0
        res <- ts(obj$residuals[-(1:GarchOrder[2])])
    }
    else if (class(obj) == "fGARCH") {
        order <- 0
        res <- ts(residuals(obj, standardize = TRUE))
    }
    else if (class(obj) == "list"){
        order <- obj$order
        if(is.null(order))
          order <- 0
        else 
          order <- order
        res <- obj$res
    }
  return(list(order = order, res = res))
}