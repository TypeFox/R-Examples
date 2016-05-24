estN.stat <-
function(z){
t <- as.integer(log(length(z) + 1) / log(2))
Mat <- matrix(0, ncol = t, nrow = 2 ^ t - 1)
for(i in 1:(2 ^ t - 1)) Mat[i,] <- rev(as.integer(intToBits(i))[1:t])
M <- sum(z)
D <- M - mean(sapply(1:t, function(j) Sub.S(z, t, Mat, j)))
Chat <- 1 - mean(sapply(1:t, function(j) Sub.S(z, t, Mat, j)) / sapply(1:t, function(j) Sub.n(z, t, Mat, j)))
    stat = cbind(M, round(D, 3), round(Chat,3))
    colnames(stat) = c("M", "D", "Chat")
    return(stat)
}
