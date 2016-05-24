mmds <- function (active, pc=3, group.file=NULL) {

 #active validation
    if (!is.matrix(active))
	stop("active is not a matrix")
    if (nrow(active) != ncol(active))
	stop("active is not square")
    if (!isSymmetric.matrix(active))
	stop("active is not symmetric")
    #if (sum(diag(active)) != 0)
    #	stop("active diagonal values are not zero")

    if (!is.numeric(active))                         
        stop("numeric values expected in active") 
    if (any(!is.finite(active)))
	stop("infinite or missing values not allowed in active")
    if (is.null(rownames(active)))
	rownames(active) <- paste("A", 1:nrow(active), sep = "")
    if (is.null(colnames(active)))
	colnames(active) <- rownames(active)

	#results will be stored in a list
	res <- list ()

	#MDS of active data
	D <- active^2

	#identity matrix
        I <- diag(1, nrow(active))
	#I <- matrix(0, nrow=nrow(active),ncol=nrow(active))

	#active matrix of ones
	ONES <- matrix(1, nrow = nrow(active), ncol = 1)
	#m vector of mass
        m<-matrix(1/nrow(active),nrow = nrow(active), ncol = 1)
	BigI<-I-(ONES%*%t(m))
	#compute active cross-product matrix
	S <-  -0.5 * (BigI %*% D %*% t(BigI))
	eigen <- eigen(S)
	res$eigen <- round(eigen$values[1:pc], 3)
	eigen.perc <- (abs(eigen$values) * 100) / sum(eigen$values[eigen$values>0])
	res$eigen.perc <- round(eigen.perc[1:pc], 3)
	#only positive eigenvalues are kept
	eigen$vectors <- eigen$vectors[, eigen$values > 0]
	eigen$values <- eigen$values[eigen$values > 0]
	
	#res$eigen <- round(eigen$values[1:pc], 3)
	res$source<-list()
	res$source$D <- D
	res$source$m<-m
	#check principal components
	if (pc < 2)
		pc <- 3
	if (pc > length(eigen$values))
		pc <- length(eigen$values)

	#eigenvalues are transformed into percentage


	#compute active matrix of factor scores
	F <- diag(as.vector(m)^(-0.5)) %*% eigen$vectors %*% diag(eigen$values^0.5)

	coord <- data.frame(F[, 1:pc])
	rownames(coord) <- rownames(active)
	colnames(coord) <- paste ("PC", (1:pc), sep = "")
	res$coord = round(coord, 3)
	res$group<-matrix(c("NoGroup","black"),1)
	colnames(res$group)<-c("group","color")
	#res$active.group<-as.data.frame.matrix(res$group)
	res$col<-matrix(c("","NoGroup","black"),1)
        colnames(res$col)<-c("element","group","color")
	#res$col<-as.data.frame.matrix(res$col)
	class (res) <- c("mmds")
	if(!is.null(group.file)){
		res<-col.group(res,group.file)
        }
	return (res)
}
