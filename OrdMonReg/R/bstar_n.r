bstar_n <- function(g1, w1, g2, w2){

n <- length(g1)

## left side of min
prod2 <- g2 * w2
c.gw2 <- rev(cumsum(rev(prod2)))
c.w2 <- rev(cumsum(rev(w2)))
tmp1 <- min(c.gw2 / c.w2)

## right side of min
tmp2 <- matrix(NA, nrow = n, ncol = n)
prod1 <- g1 * w1
c.gw1 <- rev(cumsum(rev(prod1)))
c.w1 <- rev(cumsum(rev(w1)))

for (t in 1:n){
    tbar <- t:n
    tmp2[t, tbar] <- (c.gw1[tbar] + c.gw2[t]) / (c.w1[tbar] + c.w2[t])
    } # end for
tmp2 <- min(tmp2, na.rm = TRUE)

## pooled min
res <- min(tmp1, tmp2)   

return(res)
}
