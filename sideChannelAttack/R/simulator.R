S <-
function (n) 
{
    n = intToBinary6bits(n)
    numLign = ((2^0) * n[1] + (2^1) * n[6]) + 1
    numCol = ((2^0) * n[2] + (2^1) * n[3] + (2^2) * n[4] + (2^3) * 
        n[5]) + 1
    tab = matrix(c(14, 0, 4, 15, 4, 15, 1, 12, 13, 7, 14, 8, 
        1, 4, 8, 2, 2, 14, 13, 4, 15, 2, 6, 9, 11, 13, 2, 1, 
        8, 1, 11, 7, 3, 10, 15, 5, 10, 6, 12, 11, 6, 12, 9, 3, 
        12, 11, 7, 14, 5, 9, 3, 10, 9, 5, 10, 0, 0, 3, 5, 6, 
        7, 8, 0, 13), nrow = 4)
    return(tab[numLign, numCol])
}
binaryToInt <-
function (n) 
{
    expp = length(n)
    res = 0
    j = 1
    for (i in expp:1) {
        res = res + (n[j] * 2^(i - 1))
        j = j + 1
    }
    return(res)
}
hm <-
function (n) 
{
    res = 0
    for (i in 1:length(n)) res = res + n[i]
    return(res)
}
intToBinary <-
function (n) 
{
    if (n == 0) 
        return(matrix(0, ncol = 1, nrow = 1))
    res = matrix()
    while (n >= 1) {
        ttt = n%%2
        res = cbind(matrix(ttt), res)
        n = floor(n/2)
    }
    return(res[, 1:(dim(res)[2] - 1)])
}
intToBinary6bits <-
function (n) 
{
    res = intToBinary(n)
    if (length(res) < 6) {
        for (i in ((length(res)) + 1):6) res = c(0, res)
    }
    return(res)
}
simulator.Simple1 <-
function (message, key, noise=0) 
{
	if (noise < 0) {
        stop("'noise' has to be more than 0")
        return(-1)
    }
    if (!is.matrix(message)) {
        stop("'message' has to be a matrix")
        return(-1)
    }
    if (!is.vector(key)) {
        stop("'key' has to be a vector")
        return(-1)
    }
    if ((dim(message)[2]) != 6 || (length(key)) != 6) {
        stop("the size of 'messages' and the 'key' has to be in 6 bits")
        return(-1)
    }
	if( length( which(message!=0 & message!=1)) != 0 ){
        stop("'message' has to be a vector in binary of 6 boxes. Each box has to be 1 or 0")
        return(-1)
    }
	if( length( which(key!=0 & key!=1)) != 0 ){
        stop("'key' has to be a vector in binary of 6 boxes. Each box has to be 1 or 0")
        return(-1)
    }
    tracee = matrix(ncol = (dim(message)[1]), nrow = 1)
    avantSurLeBus = rep(0, 6)
    for (t in 1:(dim(message)[1])) {
        tmp1 = binaryToInt(as.numeric(xor(message[t, ], key)))
        tracee[t] = hm(as.numeric(xor(intToBinary6bits(S(tmp1)), 
            avantSurLeBus))) + rnorm(1,mean=0,sd=noise)
        avantSurLeBus = intToBinary6bits(S(tmp1))
    }
    return(tracee)
}
