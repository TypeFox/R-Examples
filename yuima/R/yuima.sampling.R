
# general behaviour
# if grid is specified, the following are derived from it
# grid ->  n, delta, Initial, Terminal, regular, random

# grid is ALWAYS a list, possibly of dimension 1
# if it is not a list, we transform it to a list

# if grid is 1-dim, no problem, but we can have more grids. 
# in this case it is better to have a listI replace grid
# with alist

##Constructor and Initializer of class 'sampling'

# we convert objects to "zoo" internally

 
# to be fixed: the grid should always be prepared unless it is random sampling


# which.delta: check is grid is regular. If regular returns the delta, otherwise NA

which.delta <- function(x) ifelse(length(unique(round(diff(x),7)))==1, diff(x)[1], NA)

setMethod("initialize", "yuima.sampling",
function(.Object, Initial, Terminal, n, delta, grid, random, 
regular, sdelta, sgrid, oindex, interpolation){ 
				.Object@sdelta <- as.numeric(NULL) 	 
				.Object@sgrid <- as.numeric(NULL) 	 
				.Object@oindex <- as.numeric(NULL) 	 
				.Object@interpolation <- interpolation 	 
				
# grid given              
				if(!is.null(grid)){
					if(!is.list(grid)){
						yuima.warn("attempting to coerce 'grid' to a list, unexpected results may occur!")
						grid <- list(grid)
					}
				   grid <- lapply(grid, sort) # we make sure grids are ordered
				   .Object@grid <- grid
				   .Object@Initial <- sapply(grid, min)
				   .Object@Terminal <- sapply(grid, max)
				   .Object@n <- sapply(grid, function(x) length(x))
				   .Object@random <- FALSE
				   .Object@delta <- as.numeric(sapply(grid, which.delta))
				   .Object@regular <- !any(is.na( .Object@delta ) )	
				   return(.Object)
				}
# grid is missing, but random sampling					
				if(!is.logical(random)){
					.Object@regular <- FALSE
					.Object@Terminal <- Terminal	
					.Object@Initial <- Initial
					.Object@n <- numeric(0)
					.Object@delta <- numeric(0)	
					.Object@random <- random	
					return(.Object)
				}
					 


# grid is missing, but non random sampling	


                nTerm <- 0
	            if(!missing(Terminal)) nTerm <- length(Terminal)
				nInit <- 0
				if(!missing(Initial))  nInit <- length(Initial)
				nObs <- 0
				if(!missing(n)) nObs <- length(n)
				nDelta <- 0
				if(!any(is.na(delta))) nDelta <- length(delta)
				
				grid <- list()
	
# Initial + delta + n (+ Terminal: ignored) => Terminal
				if(nInit>0 & nDelta>0 & nObs>0){
				  	dims <- c(nInit, nDelta, nObs)
					ndim <- dims[ which.max(dims) ]
					Initial <- rep(Initial, ndim)[1:ndim] 	
					delta <- rep(delta, ndim)[1:ndim]
					n <- rep(n, ndim)[1:ndim] 	
					Terminal <- Initial + n*delta 	
					yuima.warn("'Terminal' (re)defined.")
					for(i in 1:ndim)
						grid[[i]] <- seq(Initial[i], Terminal[i], by=delta[i])
					
					.Object@Terminal <- Terminal
					.Object@Initial <- Initial
					.Object@n <- n	 
					.Object@delta <- delta
					.Object@grid <- grid
					.Object@random <- FALSE
					.Object@regular <- TRUE	 
					
					return(.Object)
				}

					 
# Initial + Terminal + n (+ delta: ignored) => delta
				if(nInit>0 & nTerm>0 & nObs>0){
					dims <- c(nInit, nTerm, nObs)
					ndim <- dims[ which.max(dims) ]
                    Initial <- rep(Initial, ndim)[1:ndim]
					Terminal <- rep(Terminal, ndim)[1:ndim]
					if( any(Terminal < Initial))
						stop("\nYUIMA: 'Terminal' < 'Initial'\n")
                    n <- as.integer(n)
					n <- rep(n, ndim)[1:ndim] 	
					delta <- (Terminal-Initial)/n 						 
					yuima.warn("'delta' (re)defined.")
					for(i in 1:ndim)
						grid[[i]] <- seq(Initial[i], Terminal[i], by=delta[i])					
				}

# Initial + delta + Terminal ( ignored) => n
                if(nInit>0 & nTerm>0 & nDelta>0){
                    dims <- c(nInit, nTerm, nDelta)
                    ndim <- dims[ which.max(dims) ]
                    delta <- rep(delta, ndim)[1:ndim]
                    Initial <- rep(Initial, ndim)[1:ndim]
                    Terminal <- rep(Terminal, ndim)[1:ndim]
                    if( any(Terminal < Initial))
                        stop("\nYUIMA: 'Terminal' < 'Initial'\n")
                    n <- as.integer((Terminal-Initial)/delta)
                    n <- rep(n, ndim)[1:ndim]
                    yuima.warn("'n' (re)defined.")
                    for(i in 1:ndim)
                        grid[[i]] <- seq(Initial[i], Terminal[i], by=delta[i])
                }


				.Object@Terminal <- Terminal
				.Object@Initial <- Initial
				.Object@n <- n	 
				.Object@delta <- delta
				.Object@grid <- grid
				.Object@random <- FALSE
				.Object@regular <- TRUE	 
					 
				return(.Object)
})

setSampling <- function(Initial=0, Terminal=1, n=100, delta, 
 grid, random=FALSE, sdelta=as.numeric(NULL), 
 sgrid=as.numeric(NULL), interpolation="pt" ){
  if(missing(delta))	delta <- NA
  if(missing(grid))	 grid <- NULL
  return(new("yuima.sampling", Initial=Initial, Terminal=Terminal, 
	n=n, delta=delta, grid=grid, random=random, 
			 regular=TRUE, sdelta=sdelta, sgrid=sgrid,
			 interpolation=interpolation))
}

#setMethod("initialize", "yuima.sampling",
#function(.Object, Initial, Terminal, n, delta, grid, random, 
#regular, sdelta, sgrid, oindex, interpolation){  
#	.Object@sdelta <- as.numeric(NULL) 	 
#	.Object@sgrid <- as.numeric(NULL) 	 
#	.Object@oindex <- as.numeric(NULL) 	 
#	.Object@interpolation <- interpolation 	 
## grid given                 
#	if(length(grid)>0){
#		testInitial<-(min(grid)==Initial)
#		testTerminal<-(max(grid)==Terminal)
#		testn<-(abs(n-diff(range(grid))/mean(diff(grid))+1)<10^(-10))
#		testdelta<-(abs(delta-mean(diff(grid)))<10^(-10))
#		testregular<-all(abs(diff(diff(grid)))<10^(-10))
#		
#		if(!testInitial){
#			cat("\n Start time has been set with the grid \n")
#		}
#		if(!testTerminal){
#			cat("\n Terminal time has been set with the grid \n")
#		}
#		if(!testn){
#			cat("\n Division has been set with the grid \n")
#		}
#		if(!testdelta){
#			cat("\n delta has been set with the grid \n")
#		}
#		if(testregular){
#			.Object@n <- diff(range(grid))/mean(diff(grid))+1
#			.Object@delta    <- mean(diff(grid))
#			.Object@regular  <- TRUE
#		}else{
#			.Object@n <- length(grid)-1
#			.Object@delta    <- as.numeric(NULL)
#			.Object@regular  <- FALSE  
#		}
#		.Object@grid     <- grid
#		.Object@Initial  <- min(grid)
#		.Object@Terminal <- max(grid)   
#		.Object@random   <- random   
#	}else{ 
## There is no grid
#		eqn <- length(Terminal)
#		if(length(Terminal)==length(n)){
#			.Object@Initial  <- Initial
#			.Object@Terminal <- Terminal
#			.Object@n <- n
#			.Object@delta    <- (Terminal-Initial)/n
#			.Object@grid     <- seq(Initial,Terminal,by=.Object@delta)
#			.Object@random   <- FALSE
#			.Object@regular  <- regular
#		}else if(length(Terminal)==1){               
#			.Object@Initial  <- Initial
#			.Object@Terminal <- rep(Terminal, length(n))
#			.Object@n <- n	 
#			.Object@delta    <- (Terminal-Initial)/n
#			.Object@grid     <- seq(Initial,Terminal,by=.Object@delta)
#			.Object@random   <- FALSE
#			.Object@regular  <- regular	 
#		}else if(length(n)==1){
#			.Object@Initial  <- Initial
#			.Object@Terminal <- Terminal
#			.Object@n <- rep(n, length(Terminal))
#			.Object@delta    <- (Terminal-Initial)/n
#			.Object@grid     <- seq(Initial,Terminal,by=.Object@delta)
#			.Object@random   <- FALSE
#			.Object@regular  <- regular
#		}else{
#			cat("\nDimension missmatch.\n")
#			return(NULL)
#		}}
#	return(.Object)
#})
#