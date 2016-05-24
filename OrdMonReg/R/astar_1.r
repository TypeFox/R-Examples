astar_1 <- function(g1, w1, g2, w2){

n <- length(g1)

## left side of max
prod1 <- g1 * w1
c.gw1 <- cumsum(prod1)
c.w1 <- cumsum(w1)
tmp1 <- max(c.gw1 / c.w1)

## right side of max
tmp2 <- matrix(NA, nrow = n, ncol = n)
prod2 <- g2 * w2
c.gw2 <- cumsum(prod2)
c.w2 <- cumsum(w2)

for (tbar in 1:n){tmp2[tbar:n, tbar] <- (c.gw1[tbar:n] + c.gw2[tbar]) / (c.w1[tbar:n] + c.w2[tbar])}
tmp2 <- max(tmp2, na.rm = TRUE)

## pooled min
res <- max(tmp1, tmp2) 

return(res)
}
