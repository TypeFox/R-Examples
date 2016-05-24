
performance.diagram <- function(...){

far <- seq(1,0, length = 501)### 1 - far
h   <- seq(0,1, length = 501) ### pod


f   <- function(far, h){  (1- far)*h/(1-far*(1-h) )}
g   <- function(far, h){  h/(1-far)}
hh  <- function(h, b) { h/b  }


TS <- B <- matrix(NA, nrow = 501, ncol = 501)

 
for(i in 1:501){
  for(j in 1:501){
TS[i,j] <- f(far[i], h[j])


}
}

contour(t(TS), xlim = c(0,1), ylim = c(0,1), xlab = "Success Ratio",
        ylab = "Probability of Detection", ... )
BB <- c(0.3,0.5, 0.8, 1, 1.3, 1.5,2, 3, 5, 10)


x0 <- 0
y0 <- 0
x1 <- 1
y1 <- hh(1, 1/BB)


segments(x0, y0, x1, y1, lty = 2, col = 1)

id <- y1<= 1
mtext(side = 4, text = y1[id], at = y1[id], line = 0.3, cex = 0.7, las = 2)

id <- y1> 1
mtext(side = 3, text = y1[id], at = 1/y1[id], line = 0.3, cex = 0.7) 
}  ## close function
