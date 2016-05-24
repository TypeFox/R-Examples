StereoPlane <-
function(my.strike = 000, my.dip = 45, my.color = 'black'){

    my.strike <- my.strike * (pi / 180)
    my.dip <-  my.dip * (pi / 180)
    
    my.phi <- seq(from = (-1 * pi / 2), to = (pi / 2), length = 180)
    
    my.strike.cos <- cos(my.strike)
    my.strike.sin <- sin(my.strike)
  if(my.dip != 0){
    my.lambda <- pi / 2 - my.dip
    my.alpha <- acos(cos(my.phi) * cos(my.lambda))
    my.tq <- sqrt(2) * sin(my.alpha / 2)
    my.sin.t <- sin(my.phi) / sin(my.alpha)
    my.temps <- sapply(my.sin.t, function(x){1 - x * x})
    my.temps[which(my.temps < 0)] <- 0
    my.xt <- my.tq * sqrt(my.temps)
    my.yt <- my.tq * my.sin.t
  }else{
    my.xt <- c(sin(my.phi + pi / 2), sin(my.phi - pi / 2))
    my.yt <- c(cos(my.phi + pi / 2), cos(my.phi - pi / 2))
  }
  my.x <- my.strike.cos * my.xt + my.strike.sin * my.yt
  my.y <- -1 * my.strike.sin * my.xt + my.strike.cos * my.yt
  lines(my.x, my.y, col = my.color, lwd = .5)
}
