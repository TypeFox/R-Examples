make.circle.discovery <- function(cn, diam, diam.y = NA, col = "darkgrey"){

 if (is.na(diam.y)){
  diam.y<-diam
 }
  phi <- seq(0,2 * pi,length = 1000)
  
  complex.circle <- complex(modulus = 1,argument = phi)
  
  polygon(x = cn[1] + diam * Re(complex.circle) / 2, y = cn[2] + diam.y * Im(complex.circle) / 2, border = col)
}

make.arrow.circle <- function(x, y, diam, diam.y = NA, lwd = 1, code = 2){
    
 if (is.na(diam.y)){
  diam.y <- diam
 }

  xy.factor <- c(diam,diam.y)
  x <- x / xy.factor
  y <- y / xy.factor

  diam <- 1

   if (x[1] > y[1]){
    #Calculate the skip vector v:
    alpha <- atan((x[2] - y[2])/ (x[1] - y[1]))
    v <- c(cos(alpha) * diam / 2, sin(alpha) * diam / 2)
   }
   
    if (x[1] == y[1]){
     if (x[2] > y[2]){
      v <- c(0, diam / 2)
     } else {
        v <- c(0, -diam / 2)
       }
   }

    if (x[1] <y[1]){
     #Calculate the skip vector v:
     alpha <- atan((x[2] - y[2]) / (y[1] - x[1]))
     v <- c(-cos(alpha) * diam / 2, sin(alpha) * diam / 2)
    }

    if (sum((x-y)^2)<=diam^2){
     x.new <- x + v
     y.new <- y + v
    } else {
       x.new <- x - v
       y.new <- y + v
      }

    x.new <- x.new * xy.factor
    y.new <- y.new * xy.factor
    arrows(x.new[1], x.new[2], y.new[1], y.new[2], lwd = lwd, length = 0.1, angle = 20, code = code)
}
