
##:: function subsampling
##:: takes any yuima object with data or a yuima.data object and
##:: performs subsampling according to some method

# poisson.random.sampling
# returns sample of data using poisson sampling

setGeneric("subsampling", 
function(x, sampling, ...) 
 standardGeneric("subsampling")
)

setMethod("subsampling","yuima", 
function(x, sampling, ...){
obj <- NULL
	if(missing(sampling)){
		obj <- subsampling(x@data, setSampling(...))
	} else {
		obj <- subsampling(x@data, sampling=sampling)
	}
 obj@model <- x@model
 return(obj)
}
)

setMethod("subsampling", "yuima.data",
function(x, sampling=sampling){

#tmpsamp <- NULL
#if(missing(sampling)){
#	tmpsamp <- setSampling(Initial = Initial, Terminal = Terminal, 
#				delta = delta, grid = grid, random = random, sdelta=sdelta, 
#				sgrid=sgrid, interpolation=interpolation)
#} else {
	tmpsamp <- sampling
#}


 Data <- get.zoo.data(x)
 n.data <- length(Data)
 tmpgrid <- vector(n.data, mode="list")
 tmpsamp@grid <- rep(tmpsamp@grid, n.data)[1:n.data]
	
# prepares a grid of times

	 if(is.logical(tmpsamp@random)){
      if(tmpsamp@random)
		 stop("wrong random sampling specification")
	  if(length(tmpsamp@delta) < n.data)
		 tmpsamp@delta <- rep( tmpsamp@delta, n.data)[1:n.data]
      for(i in 1:n.data){
		  tmpgrid[[i]] <- tmpsamp@grid[[i]] #seq(start(Data[[i]]), end(Data[[i]]), by=tmpsamp@delta[i])
		  tmpsamp@regular[i] <- TRUE
		  tmpsamp@random[i] <- FALSE
	  }
     }

# random sampling
	 if(is.list(tmpsamp@random)){
		 rdist <- c(tmpsamp@random$rdist)
		 if(is.null(rdist))
			stop("provide at least `rdist' argument for random sampling")
		 n.rdist <- length(rdist)
		 r.gen <- rep( rdist, n.data) # eventually reciclying arguments
		 r.gen <- r.gen[1:n.data]
		 for(i in 1:n.data){
			 tmptime <- start(Data[[i]])	
			 T <- end(Data[[i]])
			 while(	sum( tmptime ) < T )
				tmptime <- c(tmptime, r.gen[[i]](1))
			 tmpgrid[[i]] <- cumsum(tmptime)
			 if(tail(tmpgrid[[i]],1)>T)
				tmpgrid[[i]] <- tmpgrid[[i]][-length(tmpgrid[[i]])]
		 }
	 } 

# prepares original index slot	 
	 oindex <- vector(n.data, mode="list")
# checks for interpolation method, if not in the list uses "pt"
	 interpolation <- tmpsamp@interpolation
	 int.methods <- c("previous", "pt", "next", "nt", "none", "exact", 
					  "lin", "linear")
     if(! (interpolation %in% int.methods) )
	  interpolation <- "pt"

	 for(i in 1:n.data){
	  oindex[[i]] <- time(Data[[i]]) 
	  idx <- numeric(0)
	  newData <- NULL
	  lGrid <- length(tmpgrid[[i]]) 
	  if( interpolation %in% c("previous", "pt")){	 
		 idx <- as.numeric(sapply(tmpgrid[[i]], function(x) max(which(oindex[[i]] <= x))))
		 newData <- sapply(1:lGrid, function(x) as.numeric(Data[[i]][idx[x]]))  
		 oindex[[i]] <- sapply(1:lGrid, function(x) time(Data[[i]])[idx[x]])
	  }
	  if( interpolation %in% c("next", "nt")){	 
		 idx <- as.numeric(sapply(tmpgrid[[i]], function(x) min(which(oindex[[i]] >= x))))
		 newData <- sapply(1:lGrid, function(x) as.numeric(Data[[i]][idx[x]]))  
		 oindex[[i]] <- sapply(1:lGrid, function(x) time(Data[[i]])[idx[x]])
	  }
	  if( interpolation %in% c("none", "exact")){
		 idx <- match(tmpgrid[[i]], oindex[[i]])
		 newData <- sapply(1:lGrid, function(x) as.numeric(Data[[i]][idx[x]]))  
		 oindex[[i]] <- sapply(1:lGrid, function(x) time(Data[[i]])[idx[x]])
	  }

	  if( interpolation %in% c("lin", "linear")){
		 idx.l <- as.numeric(sapply(tmpgrid[[i]], function(x) max(which(oindex[[i]] <= x))))
		 idx.r <- as.numeric(sapply(tmpgrid[[i]], function(x) min(which(oindex[[i]] >= x))))
		 f.int <- function(u)
		  (as.numeric(Data[[i]][idx.r[u]])+as.numeric(Data[[i]][idx.l[u]]))/2
		 newData <- sapply(1:lGrid, f.int ) 
		 oindex[[i]] <- sapply(1:lGrid, function(u) time(Data[[i]])[idx.l[u]])
 	  }
	  Data[[i]] <- zoo(newData, order.by=tmpgrid[[i]])	 
	  tmpsamp@Terminal[i] <- end(Data[[i]])	 
	  tmpsamp@Initial[i] <- start(Data[[i]])	 
	  tmpsamp@n[i] <- length(Data[[i]])	 	 
	 }
	 

	 tmpsamp@oindex <- oindex
	 tmpsamp@grid <- tmpgrid
	 tmpsamp@regular <- sapply(1:n.data, function(x) sum(abs(diff(diff(tmpgrid[[x]]))))<1e-3)
	 tmpsamp@delta <- sapply(1:n.data, function(x) ifelse(tmpsamp@regular[x], diff(tmpgrid[[x]])[1],  numeric(0)))
	 obj <- NULL
	 tmpsamp@interpolation <- interpolation

	 x@zoo.data <- Data 		 
	 obj <- setYuima(data=x, sampling=tmpsamp)
	 return(obj)
 } ### end method
)


