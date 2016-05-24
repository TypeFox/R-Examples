`irt2ctt` <-
  ## Lord, 1980, p. 33-34  
  function (a = 1, b = 0, c = 0, d = 1, model = "LOGISTIC") 
  {
    note              <- NULL
    if (c > 0 || d < 1) 
      warning("Considering the values of c and d parameters, the rpbis value is not valid, nor the normal.parameters vector.")
    
    if (model == "LOGISTIC") { # sans la constante D
      a.irt        <- a
      a.normal     <- a*1.702    
      rpbis        <- as.numeric(a.normal/sqrt(1 + a.normal^2))   
      p            <- as.numeric(1 - pnorm(rpbis * b, mean = 0, sd = 1, lower.tail = TRUE))
    }
    
    if (model == "NORMAL") { # avec la constante D
      a.irt        <- a/1.702  
      a.normal     <- a      
      rpbis        <- as.numeric(a.normal/sqrt(1 + a.normal^2))       
      p            <- as.numeric(1 - pnorm(rpbis * b, mean = 0, sd = 1, lower.tail = TRUE))
    }
    
    parameters <- c(rpbis = rpbis, difficulty = p)
    return(parameters)
  }
