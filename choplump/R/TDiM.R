`TDiM` <-
function(S,R){ 
    N<-length(S)
    (1/sqrt(N-1))*(sum(S*R) - N*mean(S)*mean(R))/(sd(S)*sd(R))
}

