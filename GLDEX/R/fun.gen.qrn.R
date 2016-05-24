"fun.gen.qrn" <-
function(n, dimension, scrambling,FUN="runif.sobol"){

if(FUN == "runif.sobol"){
result<-runif.sobol(n = n, dimension = dimension, scrambling = scrambling)
}

if(FUN == "runif.halton"){
result<-runif.halton(n = n, dimension = dimension)
}

if(FUN == "QUnif"){
result<-QUnif(n = n, p = dimension, leap = scrambling)
}

result

}

