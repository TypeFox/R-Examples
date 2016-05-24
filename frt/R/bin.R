bin <-
function(x){
num <- 0
i <- 0
while (x > 0){
i <- i + 1
num[i] <- x %% 2
x <- floor(x / 2)
}
return(num)
}

