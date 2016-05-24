# from phangorn
fast.table <- function (data)                                                            
{                                                                                 
    if(!is.data.frame(data)) 
        data = as.data.frame(data, stringsAsFactors = FALSE)                    
    da = do.call("paste", c(data, sep = "\r"))                                    
    ind = !duplicated(da)                                                                  
    levels = da[ind]                                                                       
    cat <- factor(da,levels = levels)                                                      
    nl <- length(levels(cat))                                                        
    bin <- (as.integer(cat) - 1)                                                           
    pd <- nl                                                                               
    bin <- bin[!is.na(bin)]                                                                
    if (length(bin)) bin <- bin + 1                                                        
    y <- tabulate(bin, pd)                                                                 
    result=list(index = bin, weights = y, data = data[ind,])                                                                                  
    result                                                                                 
}                                                                                        


mydist <- function(data, k=20, distance = 2){
    m <- dim(data)[1]
    data <- cbind(data, rnorm(m, sd=1e-6))
    q <- dim(data)[2]
    k1 = k+1L
    dmtmp <- .C("dmEuclid", as.double(data), as.double(data), 
        as.integer(m), as.integer(m), as.integer(q), dm = double(m*k1), 
        cl = integer(m*k1), k = as.integer(k1), as.double(distance), 
        as.double(rep(1,q)), PACKAGE = "kknn")
    D <- matrix(dmtmp$dm, nrow = m, ncol = k1)[,-1, drop=FALSE]
    C <- matrix(dmtmp$cl + 1L, nrow = m, ncol = k1)[,-1, drop=FALSE]
    list(D, C)
}


getClosest = function(X, Y){
    m = nrow(Y)
    n = nrow(X)
    res = matrix(0, n, m)
    for(i in 1:m){
            tmp = (X-rep(Y[i,], each=n))**2
            res[,i] = rowSums(tmp)
    }
    apply(res, 2, which.min)
}


Laplacian <- function(DC, k, normalize="none"){
    normalize <- match.arg(normalize, c("none", "symmetric", "random-walk"))   
    m <- dim(DC[[1]])[1]
    INDEX = matrix(c( rep(1:m,k), as.vector(DC[[2]])) , ncol=2)
    ind = which(INDEX[,2] < INDEX[,1])
    INDEX[ind, ] = INDEX[ind, c(2,1), drop=FALSE]  
    INDEX2 = fast.table(INDEX) 
    ind = which(!duplicated(INDEX2[[1]]))
    INDEX =INDEX2[[3]]
    i = c(INDEX[,1],INDEX[,2])
    j = c(INDEX[,2],INDEX[,1])
    X = as.vector(DC[[1]])[ind]
    x =  c(X, X)
# graph.laplacian ??
    result <- sparseMatrix(i = i, j = j, x=x, dims = c(m,m))
    D = rowSums(result)
    if(normalize=="none") return(Diagonal(x=D) - result)
    if(normalize=="symmetric"){
         TMP = Diagonal(x=1/sqrt(D))
         result = TMP %*% result %*% TMP  
         return(Diagonal(m) - result)
    }
    if(normalize=="random-walk"){
        return(Diagonal(m) - Diagonal(x=1/D)%*%result)
    }
    result
}


AUC = function(y){
    l = length(y)
    x = 0:(l-1)
    y = y - y[1]
    res = numeric(0)

    for(i in 1:l){
        A = 0
        A = y[i]*(i-1)/2 
        B = y[i] * (l-i)
        C = (y[l] - y[i]) *  (l-i) / 2
        res[i] = A+B+C
    }
    res
}


specClust <- function (data, centers=NULL, nn = 7, method = "symmetric", gmax=NULL, ...) 
{
    call = match.call()
    if(is.data.frame(data)) data = as.matrix(data)
    # unique data points
    da = apply(data,1, paste, collapse="#")
    indUnique = which(!duplicated(da))
    indAll = match(da, da[indUnique])	
	
    data2 = data
    data  = data[indUnique, ]
    n <- nrow(data) 

    data = scale(data, FALSE, TRUE)


    if(is.null(gmax)){ 
         if(!is.null(centers)) gmax = centers - 1L
         else gmax = 1L
    }
    test=TRUE
    while(test){ 
        DC = mydist(data, nn)
        sif <- rbind(1:n, as.vector(DC[[2]]))
        g <- graph(sif, directed=FALSE)   
        g <- decompose(g, min.vertices=4)
        if (length(g) > 1) {
            #warning("graph not connected")
            if(length(g)>=gmax) nn = nn+2       
            else test=FALSE
        }
        else test=FALSE
    }

    W <- DC[[1]]
    n <- nrow(data) 
    wi <- W[,nn]
    SC <- matrix(1, nrow(W), nn) 
    SC[] <-  wi[DC[[2]]] * wi
    W = W^2 / SC

    alpha=1/(2*(nn+1))
    qua=abs(qnorm(alpha))
    W = W*qua
    W = dnorm(W, sd = 1)
    
    DC[[1]] = W
    L = Laplacian(DC, nn, method)

    f <- function(x, extra) as.vector(extra %*% x)

    if(is.null(centers))kmax = 25
    else kmax = max(centers)

    U <- arpack(f, extra = L, options = list(n = n, which = "SM", 
        nev = kmax, ncv = 2 * kmax, mode=1), sym = TRUE)
    ind <- order(U[[1]])
    U[[2]] = U[[2]][indAll, ind]
    U[[1]] = U[[1]][ind]
    if (is.null(centers)) {
        tmp = which.max(diff(U[[1]]))+1  
        centers = which.min(AUC(U[[1]][1:tmp]))
    }
    if(method == "symmetric"){
        rs = sqrt(rowSums(U[[2]]^2))
        U[[2]] =  U[[2]]/rs 
    }
    result = kmeans(U[[2]], centers = centers, nstart = 20, ...)
    archeType = getClosest(U[[2]][indAll, ], result$centers) 
    result$eigenvalue = U[[1]]
    result$eigenvector = U[[2]]
    result$data = data2
    result$indAll = indAll 
    result$indUnique = indUnique  
    result$L = L
    result$archetype = archeType
    result$call = call
    class(result) = c("specClust", "kmeans")
    result
}


plot.specClust = function(x, ...){
    plot(x$eigenvalue, xlab="# clusters", ylab="Eigen values", ...)
    abline(v=max(x$cluster), col="red")
}


