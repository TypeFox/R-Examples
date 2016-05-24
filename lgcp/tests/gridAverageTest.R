library(lgcp)

n <- 3 # "number of iterations"
m <- 3 # "number of aggregated time points"
temporal.fitted <- c(1,2,3) # required to set up MonteCarloAverage
SpatialOnlyMode <- FALSE

M <- 5
N <- 5

Y <- list()
for(i in 1:n){
    Y[[i]] <- list()
    for (j in 1:m){
        Y[[i]][[j]] <- matrix(runif(25),5,5)
    }
}

fun1 <- function(Y){
    return(Y)
}
fun2 <- function(Y){
    return(Y^2)
}

mca <- MonteCarloAverage(list("fun1","fun2"))

GAinitialise(mca)
oldtags <- list()
for(i in 1:n){
    oldtags$Y <- Y[[i]]
    GAupdate(mca)
}

GAfinalise(mca)

ret <- GAreturnvalue(mca)

Ymean <- as.list(rep(0,m))
Ymean2 <- as.list(rep(0,m))
for(i in 1:n){
    for (j in 1:m){
        Ymean[[j]] <- Ymean[[j]] + Y[[i]][[j]]
        Ymean2[[j]] <- Ymean2[[j]] + Y[[i]][[j]]^2
    }
}
for (j in 1:m){
    Ymean[[j]] <- Ymean[[j]]/n
    Ymean2[[j]] <- Ymean2[[j]]/n
}

for (j in 1:m){
    print(all(ret$return[[1]][[j]]==Ymean[[j]]))
    print(all(ret$return[[2]][[j]]==Ymean2[[j]]))
}

