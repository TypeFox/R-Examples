DefineGame <-
function(n,V) {

	# Construyo matriz binaria
	
	vec <- c(0, 1)
	lst <- lapply(numeric(n), function(x) vec)
	Amat <- as.matrix(expand.grid(lst))[-1, ]
    BinMat <- Amat
	C <- ncol(Amat)
	F <- nrow(Amat)

	# recodifico matriz
	
	for (i in seq(1, F)) {
		for (j in seq(1, C)) {
			if (Amat[i, j] == 1) 
				Amat[i, j] <- j

		}
	}

	###################
	matrix.text <- function(mtext, sep = " ", collapse = NULL) {
		if (is.null(collapse)) 
			apply(mtext, 1, paste, collapse = sep)
		else paste(apply(mtext, 1, paste, collapse = sep), 
			collapse = collapse)
	}
	Hmat <- matrix.text(Amat)

	# quito los 0's
	
	Hmat <- Hmat <- gsub(0, "", Hmat)

	# quito los espacios
	
	Hmat <- Hmat <- gsub(" ", "", Hmat)

	# convierto a numerico y mantengo la matriz
	
	Bmat <- as.matrix(as.numeric(Hmat))

	# aplico orden lexicografico
	Lmat <- as.matrix(sort(Bmat))
    
    CoB<-data.frame(rep(0,2^n-1),row.names=Bmat)
    colnames(CoB)=c("v(i)")
    
    CoL<-data.frame(V,row.names=Lmat)
    colnames(CoL)=c("v(i)")
    
    Nam<-rownames(CoL)
    d<-length(Nam)
    for (i in seq(1,d)){
    	CoB[Nam[i],]<-CoL[Nam[i],]
    	
    }
    
    Output<-list(Binary=CoB,Lex=CoL,Ag=n,BinMat=BinMat)
    class(Output)<-"Game"
	return(Output)
}
