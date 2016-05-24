accuracy <-
function(x, y, z){
  
  acc <- sqrt(x^2 + y^2 + z^2)
  var <- mean((1-acc)^2)
  std <- sqrt(var)
  
  accuracy <- list("acc" = acc, "var" = var, "std" = std)
  
  return(accuracy)
}

accuracy.improve <- 
  function(x, y, z){
    
    x1 <- sqrt(1-x^2)
    x2 <- 0
    x3 <- x
    
    y1 <- (y/sqrt(1-x^2)) * (-x)
    y2 <- sqrt((1-x^2-y^2)/(1-x^2))
    y3 <- (y/sqrt(1-x^2)) * (sqrt(1-x^2))
    
    z1 <- y2 * (-x)
    z2 <- -(y/sqrt(1-x^2))
    z3 <- y2 * x1
    
    z.est <- z3
    z.resid <- z3-z
    acc.input <- accuracy(x,y,z)
    acc.improve <- accuracy(x,y,z3) 
        
    acclist <- list("z.est" = z.est, "z.resid" = z.resid, "acc.input" = acc.input$acc, "acc.input.var" = acc.input$var,
                    "acc.input.std" = acc.input$std, "acc.improve" = acc.improve$acc, "acc.improve.var" = acc.improve$var,
                    "acc.improve.std" = acc.improve$std)
    
    return(acclist)
    
  }