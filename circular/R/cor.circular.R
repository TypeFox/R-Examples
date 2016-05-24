
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

###############################################################
#                                                             #
#       R port: Claudio Agostinelli  <claudio@unive.it>       #
#                                                             #
#       Date: May, 26, 2006                                   #
#       Version: 0.3-1                                        #
#                                                             #
###############################################################

cor.circular <- function(x, y=NULL, test = FALSE) {

   if (!is.null(y) & NROW(x)!=NROW(y)) 
       stop("x and y must have the same number of observations")
   if (is.null(y) & NCOL(x)<2)
       stop("supply both x and y or a matrix-like x")
 
   ncx <- NCOL(x)
   ncy <- NCOL(y)
        
   # Handling missing values
   if (is.null(y)) {
       ok <- complete.cases(x)
       x <- x[ok,]
   } else {
       ok <- complete.cases(x, y)
       if (ncx==1) {
           x <- x[ok]
       } else {
           x <- x[ok,]
       }
       if (ncy==1) {
           y <- y[ok]
       } else {
           y <- y[ok,]
       }
   }
         
   n <- NROW(x)
   if (n==0) {
       warning("No observations (at least after removing missing values)")
       return(NULL)
   }
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
   if (!is.null(y)) { 
       y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
       attr(y, "class") <- attr(y, "circularp") <-  NULL
   }
        
   if (is.null(y)) {      
       result <- matrix(1, ncol=ncx, nrow=ncx)
       if (test) {
           test.stat <- matrix(0, ncol=ncx, nrow=ncx)
           p.value <- matrix(0, ncol=ncx, nrow=ncx)
       }
       for (i in 1:ncx) {
            for (j in i:ncx) {
                 res <- CorCircularRad(x=x[,i], y=x[,j], test=test)
                 result[i,j] <- result[j,i] <- res[1]
                 if (test) {
                     if (i==j) {
                         test.stat[i,i] <- NA
                         p.value[i,i] <- NA
                     } else { 
                         test.stat[i,j] <- test.stat[j,i] <- res[2] 
                         p.value[i,j] <- p.value[j,i] <- res[3]
                     }
                 }
            }
        }            
   } else {
       attributes(x) <- c(attributes(x), list(dim=c(n, ncx)))
       attributes(y) <- c(attributes(y), list(dim=c(n, ncy)))
       result <- matrix(1, ncol=ncy, nrow=ncx)
       if (test) {
           test.stat <- matrix(0, ncol=ncy, nrow=ncx)
           p.value <- matrix(0, ncol=ncy, nrow=ncx)
       }
       for (i in 1:ncx) {
            for (j in 1:ncy) {
                 res <- CorCircularRad(x=x[,i], y=y[,j], test=test)
                 result[i,j] <- res[1]
                 if (test) {
                     test.stat[i,j] <- res[2] 
                     p.value[i,j] <- res[3]
                 }
            }
       }
   }
   if (ncx==1 | (!is.null(y) & ncy==1)) {
       result <- c(result)
       if (test) {
           test.stat <- c(test.stat)
           p.value <- c(p.value)
       }
   }
   
   if (test) {
       result <- list(cor=result, statistic=test.stat, p.value=p.value)
   }
   return(result)
}

CorCircularRad <- function(x, y, test=FALSE) {
   n <- length(x)
   x.bar <- MeanCircularRad(x)
   y.bar <- MeanCircularRad(y)
   num <- sum(sin(x - x.bar) * sin(y - y.bar))
   den <- sqrt(sum(sin(x - x.bar)^2) * sum(sin(y - y.bar)^2))
   result <- num/den
   if (test) {
       l20 <- mean.default(sin(x - x.bar)^2)
       l02 <- mean.default(sin(y - y.bar)^2)
       l22 <- mean.default((sin(x - x.bar)^2) * (sin(y - y.bar)^2))
       test.stat <- sqrt((n * l20 * l02)/l22) * result
       p.value <- 2 * (1 - pnorm(abs(test.stat)))
       result <- c(result, test.stat, p.value)
   }
   return(result)
}


