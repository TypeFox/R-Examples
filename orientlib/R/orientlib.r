require(methods)
if (R.Version()$minor >= 9.0 || R.Version()$major > 1) require(stats)

# The orientation class is an abstract class; the different representations
# below descend from it

setClass('orientation')
setIs('orientation', 'vector')

setMethod('%*%', c(x = 'orientation', y = 'orientation'),
    function(x, y) {
	x <- as(x, 'rotmatrix')@x
	lx <- dim(x)[3]
	y <- as(y, 'rotmatrix')@x
	ly <- dim(y)[3]
	if (max(lx,ly) %% min(lx,ly) != 0) 
	    warning('longer object length is not a multiple of shorter object length')
	result <- array(NA, c(3,3,max(lx,ly)))
	for (i in 1:max(lx,ly)) {
	    xi <- (i-1) %% lx + 1
	    yi <- (i-1) %% ly + 1
	    result[,,i] <- x[,,xi] %*% y[,,yi]
	}
	new('rotmatrix', x = result)
    }
)

setMethod('^', c(e1 = 'orientation', e2 = 'numeric'),
    function(e1, e2) {
	x <- as(e1, 'skewvector')@x
	lx <- nrow(x)
	y <- e2
	ly <- length(y)
	l <- max(lx, ly)
	if ( l %% min(lx, ly) != 0)
	    warning('longer object length is not a multiple of shorter object length')
	result <- matrix(NA, max(lx, ly), 3)
	xi <- (0:(l-1)) %% lx + 1
	yi <- (0:(l-1)) %% ly + 1
	new('skewvector', x = x[xi, ,drop=FALSE]*y[yi])
    }
)

setMethod('t', 'orientation',
    function(x) x^(-1)
)

setMethod('mean', 'orientation',
    function(x) {
	x <- rotmatrix(x)@x
	nearest.SO3(apply(x,c(1,2),mean))
    }
)

setMethod('weighted.mean', c(x = 'orientation', w = 'numeric'),
    function(x,w) {
	stopifnot(length(x) == length(w))
	x <- rotmatrix(x)@x
	w <- array(rep(w, rep(9,length(w))),c(3,3,length(w)))
	nearest.SO3(apply(x*w,c(1,2),sum))
    }
)

setAs('matrix', 'orientation',
    def = function(from) rotmatrix(from) )
    
setAs('array', 'orientation',
    def = function(from) rotmatrix(from) ) 
    
# The rotmatrix class is an 3 x 3 x n array holding rotation matrices

setClass('rotmatrix', representation(x = 'array'))
setIs('rotmatrix','orientation')

rotmatrix <- function(a) {
    d <- dim(a)
    if (length(d) < 3) d <- c(d,1)
    a <- array(a, d)
    stopifnot(dim(a)[1] == 3, dim(a)[2] == 3)
    new('rotmatrix', x = a)
}

setMethod('rotmatrix', 'orientation',
    function(a) as(a, 'rotmatrix')
)   
    
setMethod('[', 'rotmatrix',
    function(x, i, j, ..., drop) rotmatrix(x@x[,,i,drop=FALSE])
)

setMethod('[<-', 'rotmatrix',
    function(x, i, j, value) {
	x <- x@x
	v <- as(value, 'rotmatrix')@x
	x[,,i] <- v
	new('rotmatrix', x = x)
    }
)

setMethod('[[', 'rotmatrix',
    def = function(x, i, j, ...) x@x[,,i]
)

setMethod('[[<-', 'rotmatrix',
    def = function(x, i, j, value) {
	x@x[,,i] <- value
	x
    }
)

setMethod('length', 'rotmatrix', function(x) dim(x@x)[3])

# The rotvector class is an n x 9 matrix holding the vectorized rotation
# matrices

setClass('rotvector', representation(x = 'matrix'))
setIs('rotvector', 'orientation')

rotvector <- function(m) {
    if (is.null(dim(m))) m <- matrix(m,1)
    stopifnot(ncol(m) == 9)
    new('rotvector', x = m)
}

setMethod('rotvector', 'orientation',
    function(m) as(m, 'rotvector')
)

setAs('orientation', 'rotvector',
    def = function(from) {
	x <- as(from, 'rotmatrix')
	new('rotvector', x = cbind(x@x[1,1,],x@x[2,1,],x@x[3,1,],
				   x@x[1,2,],x@x[2,2,],x@x[3,2,],
				   x@x[1,3,],x@x[2,3,],x@x[3,3,]))
    }
)

setAs('rotvector', 'rotmatrix',
    def = function(from) rotmatrix(array(t(from@x), c(3,3,nrow(from@x)))))
    
setMethod('[', 'rotvector',
    function(x, i, j, ..., drop) rotvector(x@x[i,,drop=FALSE])
)   

setMethod('[<-', 'rotvector',
    function(x, i, j, value) {
	x <- x@x
	v <- as(value, 'rotmatrix')@x
	x[i,] <- v
	new('rotvector', x = x)
    }
)

setMethod('[[', 'rotvector',
    def = function(x, i, j, ...) x@x[i,]
)

setMethod('[[<-', 'rotvector',
    def = function(x, i, j, value) {
	x@x[i,] <- value
	x
    }
)

setMethod('length', 'rotvector', function(x) nrow(x@x))

# The eulerzyx class is an n x 3 array holding Euler angles for rotations in the order
# Z, Y, X

setClass('eulerzyx', representation(x = 'matrix'))
setIs('eulerzyx', 'orientation')

eulerzyx <- function(psi, theta, phi) {
    if (missing(theta) && missing(phi)) x <- psi
    else x <- cbind(psi, theta, phi)
    if (is.null(dim(x))) x <- matrix(x,1)
    stopifnot(ncol(x) == 3)
    colnames(x) <- c('psi', 'theta', 'phi')
    new('eulerzyx', x = x)
}

setMethod('eulerzyx', c('orientation','missing','missing'),
    function(psi) as(psi, 'eulerzyx')
)

setAs('matrix', 'eulerzyx',
    def = function(from, to) {
    	stopifnot(ncol(from) == 3)
    	new('eulerzyx', x = from)
    }
)

setAs('orientation', 'eulerzyx',
    def = function(from) {
	x <- as(from, 'rotmatrix')@x
    	C <- 1-x[1,3,]^2
    	zeros <- zapsmall(C) == 0
    	theta <- ifelse(!zeros, asin(-x[1,3,]), pi/2)
    	psi <- ifelse(!zeros, atan2(x[1,2,]/C, x[1,1,]/C), 0)
    	phi <- ifelse(!zeros, atan2(x[2,3,]/C, x[3,3,]/C), 
    			 ifelse(x[1,3,] < 0, atan2(x[2,1,], x[3,1,]),
    			 		     atan2(-x[2,1,],-x[3,1,]))) 			 		    
	eulerzyx(psi, theta, phi)
    }
)

setAs('eulerzyx', 'rotmatrix',
    def = function(from) {
	x <- from@x
    	psi <- x[,1]    
    	theta <- x[,2]    
    	phi <- x[,3]
	rotmatrix(aperm(
	    array(c(cos(psi)*cos(theta),  
                    -sin(psi)*cos(phi) + cos(psi)*sin(theta)*sin(phi),
                    sin(psi)*sin(phi) + cos(psi)*sin(theta)*cos(phi),  
                    	sin(psi)*cos(theta),  
                    	cos(psi)*cos(phi) + sin(psi)*sin(theta)*sin(phi),
                    	-cos(psi)*sin(phi) + sin(psi)*sin(theta)*cos(phi),  
                    	    -sin(theta),
                    	    cos(theta)*sin(phi),  
                    	    cos(theta)*cos(phi)), c(nrow(x),3,3)),
                    c(2,3,1)))
    }
)                    

setMethod('[', 'eulerzyx',
    function(x, i, j, ..., drop) eulerzyx(x@x[i,,drop=FALSE])
)    

setMethod('[<-', 'eulerzyx',
    function(x, i, j, value) {
	x <- x@x
	v <- as(value, 'eulerzyx')@x
	x[i,] <- v
	new('eulerzyx', x = x)
    }
)

setMethod('[[', 'eulerzyx',
    def = function(x, i, j, ...) x@x[i,]
)

setMethod('[[<-', 'eulerzyx',
    def = function(x, i, j, value) {
	x@x[i,] <- value
	x
    }
)

setMethod('length', 'eulerzyx', function(x) nrow(x@x))

# The quaternion class is an n x 4 array holding unit quaternions

setClass('quaternion', representation(x = 'matrix'))
setIs('quaternion', 'orientation')

quaternion <- function(m) {
    if (is.null(dim(m))) m <- matrix(m, 1)
    stopifnot(ncol(m) == 4)
    new('quaternion', x = m)
}

setMethod('quaternion', 'orientation',
    function(m) as(m, 'quaternion')
)

setAs('quaternion', 'rotmatrix',
    def = function(from) {
	x <- from@x
     	rotmatrix(aperm(array(c(1-2*x[,2]^2-2*x[,3]^2, 
    		      2*x[,1]*x[,2]-2*x[,3]*x[,4],
    		      2*x[,1]*x[,3]+2*x[,2]*x[,4],
    		      		2*x[,1]*x[,2]+2*x[,3]*x[,4],
    		      		1-2*x[,1]^2-2*x[,3]^2,
    		      		2*x[,2]*x[,3] - 2*x[,4]*x[,1],
    		      			2*x[,1]*x[,3] - 2*x[,2]*x[,4],
    		      			2*x[,2]*x[,3] + 2*x[,1]*x[,4],
    		      			1 - 2*x[,1]^2 - 2*x[,2]^2), c(nrow(x),3,3)),
    		      c(2,3,1)))
    		      
     }
)	

setAs('orientation', 'quaternion',
    def = function(from) {
	nicesqrt <- function(x) sqrt(pmax(x,0))
	x <- as(from, 'rotmatrix')@x
    	q4 <- nicesqrt((1 + x[1,1,] + x[2,2,] + x[3,3,])/4)  # may go negative by rounding
    	zeros <- zapsmall(q4) == 0 
    	q1 <- ifelse(!zeros, (x[2,3,] - x[3,2,])/4/q4, nicesqrt(-(x[2,2,]+x[3,3,])/2))
    	q2 <- ifelse(!zeros, (x[3,1,] - x[1,3,])/4/q4, 
    				 ifelse(zapsmall(q1) != 0, x[1,2,]/2/q1, nicesqrt((1-x[3,3,])/2)))
    	q3 <- ifelse(!zeros, (x[1,2,] - x[2,1,])/4/q4, 
    				 ifelse(zapsmall(q1) != 0, x[1,3,]/2/q1, 
    				 	ifelse(zapsmall(q2) != 0, x[2,3,]/2/q2, 1)))
    	quaternion(cbind(q1,q2,q3,q4))
    }
)

setMethod('[', 'quaternion',
    function(x, i, j, ..., drop) quaternion(x@x[i,,drop=FALSE])
)    

setMethod('[<-', 'quaternion',
    function(x, i, j, value) {
	x <- x@x
	v <- as(value, 'quaternion')@x
	x[i,] <- v
	new('quaternion', x = x)
    }
)

setMethod('[[', 'quaternion',
    def = function(x, i, j, ...) x@x[i,]
)

setMethod('[[<-', 'quaternion',
    def = function(x, i, j, value) {
	x@x[i,] <- value
	x
    }
)

setMethod('length', 'quaternion', function(x) nrow(x@x))

# The skewvector class is an n x 3 matrix holding the elements of the 
# skew symmetric matrix representation of an orientation

setClass('skewvector', representation(x = 'matrix'))
setIs('skewvector', 'orientation')

skewvector <- function(m) {
    if (is.null(dim(m))) m <- matrix(m, 1)
    stopifnot(ncol(m) == 3)
    new('skewvector', x = m)
}

setMethod('skewvector', 'orientation',
    function(m) as(m, 'skewvector')
)

setAs('skewvector', 'quaternion',
    def = function(from) {
	x <- from@x
        theta <- apply(x, 1, function(row) sqrt(sum(row^2)))
        w <- cos(theta/2)
        x <- t(apply(cbind(x,w), 1, function(row) {
	    k <- max(1-row[4]^2,0)
	    sumsq <- sum(row[1:3]^2)
	    row[1:3] <- row[1:3]*ifelse(sumsq > 1.e-8, sqrt(k/sumsq), 0.5)
	    row
	}))
	quaternion(x)
    }
)

setAs('skewvector', 'rotmatrix',
    def = function(from) {
	as(as(from, 'quaternion'), 'rotmatrix')    		      
     }
)	

setAs('orientation', 'skewvector',
    def = function(from) {
	x <- as(from, 'quaternion')@x
    	theta <- 2*acos(x[,4])
    	skewvector(x[,-4]*ifelse(abs(theta)>1.e-4, theta/sqrt(pmax(1-x[,4]^2,0)), 1))
    }
)

setMethod('[', 'skewvector',
    function(x, i, j, ..., drop) skewvector(x@x[i,,drop=FALSE])
)    

setMethod('[<-', 'skewvector',
    function(x, i, j, value) {
	x <- x@x
	v <- as(value, 'skewvector')@x
	x[i,] <- v
	new('skewvector', x = x)
    }
)

setMethod('[[', 'skewvector',
    def = function(x, i, j, ...) x@x[i,]
)

setMethod('[[<-', 'skewvector',
    def = function(x, i, j, value) {
	x@x[i,] <- value
	x
    }
)

setMethod('length', 'skewvector', 
	function(x) nrow(x@x)
)

# The skewmatrix class is a 3 x 3 x n matrix holding the 
# skew symmetric matrix representation of an orientation

setClass('skewmatrix', representation(x = 'array'))
setIs('skewmatrix', 'orientation')

skewmatrix <- function(a) {
    d <- dim(a)
    if (length(d) < 3) d <- c(d,1)
    a <- array(a, d)
    stopifnot(dim(a)[1] == 3, dim(a)[2] == 3)
    new('skewmatrix', x = a)
}

setMethod('skewmatrix', 'orientation',
    function(a) as(a, 'skewmatrix')
)

setAs('skewmatrix', 'skewvector',
    def = function(from) {
	x <- from@x
	x <- cbind(x[2,3,],x[3,1,],x[1,2,])
	skewvector(x)
    }
)

setAs('skewmatrix', 'rotmatrix',
    def = function(from) {
	as(as(from, 'skewvector'), 'rotmatrix')    		      
     }
)	

setAs('orientation', 'skewmatrix',
    def = function(from) {
	x <- skewvector(from)@x
	n <- nrow(x)
    	x <- aperm(array(c(rep(0,n), x[,3], -x[,2], -x[,3], rep(0,n), x[,1], x[,2], -x[,1], rep(0,n)), c(n, 3, 3)))
    	skewmatrix(x)
    }
)

setMethod('[', 'skewmatrix',
    function(x, i, j, ..., drop) skewmatrix(x@x[,,i,drop=FALSE])
)    

setMethod('[<-', 'skewmatrix',
    function(x, i, j, value) {
	x <- x@x
	v <- skewmatrix(value)@x
	x[,,i] <- v
	skewmatrix(x)
    }
)

setMethod('[[', 'skewmatrix',
    def = function(x, i, j, ...) x@x[,,i]
)

setMethod('[[<-', 'skewmatrix',
    def = function(x, i, j, value) {
	x@x[,,i] <- value
	x
    }
)

setMethod('length', 'skewmatrix', function(x) dim(x@x)[3])


# The eulerzxz class is an n x 3 array holding Euler angles for rotations in the order
# Z, X, Z (the "x-convention")

setClass('eulerzxz', representation(x = 'matrix'))
setIs('eulerzxz', 'orientation')

eulerzxz <- function(phi, theta, psi) {
    if (missing(theta) && missing(psi)) x <- phi
    else x <- cbind(phi, theta, psi)
    if (is.null(dim(x))) x <- matrix(x,1)
    stopifnot(ncol(x) == 3)
    colnames(x) <- c('phi', 'theta', 'psi')
    new('eulerzxz', x = x)
}

setMethod('eulerzxz', c('orientation','missing','missing'),
    function(phi) as(phi, 'eulerzxz')
)

setAs('matrix', 'eulerzxz',
    def = function(from) {
    	stopifnot(ncol(from) == 3)
    	new('eulerzxz', x = from)
    }
)

setAs('orientation', 'eulerzxz',
    def = function(from) {
	x <- as(from, 'rotmatrix')@x
    	C <- x[3,2,]^2 + x[3,1,]^2
    	zeros <- zapsmall(C) == 0
    	phi <- ifelse(!zeros, atan2(x[3,1,], -x[3,2,]), atan2(x[1,2,],x[1,1,]))
    	theta <- ifelse(!zeros, atan2(ifelse(zapsmall(sin(phi)) != 0, 
    										 x[3,1,]/sin(phi), 
    										 -x[3,2,]/cos(phi)), 
    							      x[3,3,]),
    							ifelse(x[3,3,]>0, 0, pi))
    	psi <- ifelse(!zeros, atan2(x[1,3,], x[2,3,]), 0) 
	eulerzxz(phi, theta, psi)
    }
)

setAs('eulerzxz', 'rotmatrix',
    def = function(from) {
	x <- from@x
    	phi <- x[,1]    
    	theta <- x[,2]    
    	psi <- x[,3]
	rotmatrix(aperm(
	    array(c(cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),  
                -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                sin(theta)*sin(phi),  
                    	cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),  
                    	-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    	-sin(theta)*cos(phi),
                    	    sin(psi)*sin(theta),
                    	    cos(psi)*sin(theta),  
                    	    cos(theta)), c(nrow(x),3,3)),
                    c(2,3,1)))
    }
)                    

setMethod('[', 'eulerzxz',
    function(x, i, j, ..., drop) eulerzxz(x@x[i,,drop=FALSE])
)    

setMethod('[<-', 'eulerzxz',
    function(x, i, j, value) {
	x <- x@x
	v <- as(value, 'eulerzxz')@x
	x[i,] <- v
	new('eulerzxz', x = x)
    }
)

setMethod('[[', 'eulerzxz',
    def = function(x, i, j, ...) x@x[i,]
)

setMethod('[[<-', 'eulerzxz',
    def = function(x, i, j, value) {
	x@x[i,] <- value
	x
    }
)

setMethod('length', 'eulerzxz', function(x) nrow(x@x))

