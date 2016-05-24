concatList = function (lis, sep=""){
    out=lis[[1]]
    i=2
    while (i<=length(lis)){
        out=out%+%sep%+%lis[[i]]
        i=i+1
    }
    out
}



getMidPoints=function(x){
    ((c(0,x)+c(x,0))/2)[2:length(x)] 
}
#getMidPoints(1:10)


# the input is a list, typically the output from a sapply call that should be matrix, but have different length
fill.jagged.array=function(a) {
    # don't check is.matrix, because for some reason, it will return true
    max.len=max(sapply(a, length))
    sapply(a, function (e) {
        c(e, rep(NA, max.len-length(e)))
    })    
}


last = function (x, n=1, ...) {
    if (length(x)==1 & is.character(x)) tail (readLines(x), n=n, ...) # read file, length(x)==1 is needed b/c otherwise last(c("a","b")) won't behave properly
    else if (is.vector(x)) x[length(x)]
    else if (is.array(x)) x[length(x)]
    else if (is.list(x)) x[[length(x)]]
    else stop ("last(): x not supported")
}

# return a subset of data that is 1 row every thin.factor rows
ThinRows = function (dat, thin.factor=10) {
    NumRows = nrow(dat)
    dat[1:(NumRows/thin.factor)*thin.factor,]
}
thin.rows=ThinRows

#mix two arrays in an interlacing way
mix = function (a, b) {
    if (length(a)!=length(b)) print ("Length of two arguments to mix function not equal.")
    out = rep (a, each=2)
    for (i in 1:length(a)) {
        out[2*i]=b[i]
    }
    out
}


# like lag, move vector to the right/left by given number of steps
# x is a vector
shift.right = function (x, k=1) {
    p=length(x)
    x[(1+k):p]=x[1:(p-k)]
    x[1:k]=NA
    x
}

shift.left = function (x, k=1) {
    p=length(x)
    x[1:(p-k)] = x[(1+k):p]
    x[(1+p-k):p]=NA
    x
}



# trace
tr=function(m) sum(diag(m))


# serial covariance matrix
AR1 = function (p, w) {
    m = matrix(1, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            m [i,j]=w**abs(i-j)
        }
    }
    m
}

# exchangeable covariance matrix
EXCH = function (p, rho) {
    m = matrix(1, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            if (i!=j) m [i,j]=rho
        }
    }
    m
}

getUpperRight = function (matri, func=NULL) {
    n=nrow (matri)
    out= numeric ( (n-1)*n/2 )
    index=0
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            index=index+1
            out[index]=matri[i,j]
        }
    }
    if (is.null(func)) {
        out
    } else {
        func(out)
    }
}


#repeat a matrix in a block diagonal fashion
rep.matrix.block = function (x, times=2, ...) {
    orig.dim = nrow (x)     
    m = matrix (0, orig.dim * times, orig.dim * times)
    for (i in 1: times) {
        m[(1+(i-1)*orig.dim):(i*orig.dim), (1+(i-1)*orig.dim):(i*orig.dim)]=x
    }
    m
}


#it does not work on data.frame
rep.matrix = function (x, times=1, each=1, by.row=TRUE, ...) {
    if (times==0) return(NULL)
    if (by.row) {
        colnames.=colnames(x)
        new.matrix=matrix(0, nrow(x)*each*times, ncol(x) )
        for (i in 1:nrow(x)) {
            for (j in 1:each) {
                new.matrix[(i-1)*each + j,] = x[i,]
            }
        }
        
        if(times>1) {
            for (i in 2:times) {
                new.matrix[ ((i-1)*nrow(x)*each+1) : (i*nrow(x)*each), ] = new.matrix[1:(nrow(x)*each), ]
            }
        }
    
        dimnames(new.matrix)[[2]]=colnames.
        new.matrix
    }
    else {
        t ( rep.matrix(t(x), times, each, by.row=TRUE) )
    }    
}

# rep.data.frame(chi[21,], 2)
rep.data.frame = function (x, times=1, ...){
    out = x
    if (time==1) return (out)
    for (i in 2:times) {
        out = rbind (out, x)
    }
    out
}
