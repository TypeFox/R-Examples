`ctt2irt` <-
  ## Lord, 1980, p. 33-34  
  function (rpbis = 0.7071068, difficulty = 0.5) 
  {
    note <- noquote("For the moment, c and d parameters don't seem possible to be recovered from p and rpbis. These models cannot be compared for the moment.")
    t <- qnorm(1 - difficulty, mean = 0, sd = 1, lower.tail = TRUE)
    a <- rpbis/sqrt(1 - rpbis^2)
    b <- t/rpbis
    parameters.normal <- c(a = a, b = b)
    parameters.irt    <- c(a = (a / 1.702), b = b)
    result <- list(note = note, normal.parameters = parameters.normal, 
                   irt.parameters = parameters.irt)
    return(result)
  }
