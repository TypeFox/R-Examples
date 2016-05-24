gini <-
function (x, norm = FALSE, lc = FALSE, lcx = "% of objects", lcy = "% of regarded variable", lctitle = "Lorenz curve", lcg = FALSE, lcgn = FALSE) {   

x_sort <- sort(x)   
i <- length(x)   
sum_x <- sum(x_sort)
a_i <- x_sort/sum(x)
y_i <- cumsum(a_i)
z_i <- 0

for (j in 1:i) {   
z_i[j] <- j/i
}

if (lc == TRUE) {
k <- c(0,1)
l <- c(0,1)
start_start <- c(0.0, z_i[1])
start_end <- c(0.0, y_i[1])
plot(k,l, type="l", col="blue",xlab=lcx, ylab=lcy, main=lctitle)
lines(z_i, y_i)
lines (start_start, start_end)
} else

j <- 0
sum_y_i <- 0
sum_y_i[1] <- y_i[1]+0

for (j in 2:i) {
sum_y_i[j] <- y_i[j-1]+y_i[j]
}

G <- 1-1/i*sum(sum_y_i)   
G.norm <- (i/(i-1))*G   

if (lc == TRUE) {
if (lcg == TRUE) {
text (0,1,paste("G =", G), pos=4)
}
if (lcgn == TRUE) {
text (0,0.95,paste("G* =", G.norm),pos=4)
}
}

if (norm == FALSE) {   
return(G)   
} 
else 
{
return (G.norm)   
}
}
