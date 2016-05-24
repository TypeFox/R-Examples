.my.det <- function(n, index){

    return((2 * n - 3) ^ 2 + 8 * (n - index));
}

.raw.i <- function(n, index){

    return(n + 1/2 - sqrt(.my.det(n,index)) / 2);
}

.calc.j <- function(i,n,index){

    return(index + n + i * (0.5 * (i + 1) - n));
}

diss.index <- function(index, N){

    i <- floor(.raw.i(N, index));

    return(c(i = i,
             j = .calc.j(i,N,index)));
}
