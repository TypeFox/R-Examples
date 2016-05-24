
	
	
	
sortandcut <- function(x,iter=20, tau0 = NULL, fun = "BCC", method = "WBCI"){
	
	n <- nrow(x)
	m <- ncol(x)
	rlist <- list(1:n)
	clist <- list(1:m)
	
	ncl <- 1
	k <- 0
	
	if(storage.mode(x) %in% c("numeric","double")){
		ncx <- nchar(x%%1)
		if(any(ncx > 1)){
			k <- min(4,max(ncx)-1)
			xlist <- list( round( x * 10^k, 0) )
		}else{
			storage.mode(x) <- "integer"
			xlist = list(x)
		}
	}else{
		xlist = list(x)
	}
	
	
	while(k < ncl){
		if(!min(dim(xlist[[k+1]])) == 1){
			
					
		if( fun %in% c("optile", "BCC") ){
			ox <- optile(xlist[[k+1]],iter=iter)
		}
		if( fun %in% c("bary","barysort") ){
			#ox <- barysort(xlist[[k+1]])
			ox <- optile(xlist[[k+1]],iter=iter, fun = "barysort", foreign = ".Call")
		}
		if( fun %in% c("pre", "preclass") ){
			ox <- optile(xlist[[k+1]],iter=iter, fun = "preclass", foreign = ".Call")
		}
		
		rlist[[k+1]] <- rlist[[k+1]][attr(ox,"orders")[[1]]]
		clist[[k+1]] <- clist[[k+1]][attr(ox,"orders")[[2]]]

		cc <- cfluctile(ox,maxsplit=1,tau0 = tau0, plot = FALSE,method=method)
		ss <- length(attr(cc,"orders")[[1]]) > 1
		
		}else{
			# matrix has only one row and/or column
			ss <- FALSE
		}
		if(ss){
			#new split, k remains unchanged
			ncl <- ncl+1
			ox1 <- ox[attr(cc,"orders")[[1]][[1]],attr(cc,"orders")[[2]][[1]],drop=FALSE]
			ox2 <- ox[attr(cc,"orders")[[1]][[2]],attr(cc,"orders")[[2]][[2]],drop=FALSE]
	
			rs1 <- rlist[[k+1]][attr(cc,"orders")[[1]][[1]]]
			rs2 <- rlist[[k+1]][attr(cc,"orders")[[1]][[2]]]
			cs1 <- clist[[k+1]][attr(cc,"orders")[[2]][[1]]]
			cs2 <- clist[[k+1]][attr(cc,"orders")[[2]][[2]]]
			
			if(k == 0){
				xlist <- c(list(ox1),list(ox2),xlist[-1])
				rlist <- c(list(rs1),list(rs2),rlist[-1])
				clist <- c(list(cs1),list(cs2),clist[-1])
			}else{
				xlist <- c(xlist[1:k],list(ox1),list(ox2),xlist[-c(1:(k+1))])
				rlist <- c(rlist[1:k],list(rs1),list(rs2),rlist[-c(1:(k+1))])
				clist <- c(clist[1:k],list(cs1),list(cs2),clist[-c(1:(k+1))])
				
			}
			


		}else{
			k <- k+1
			
		}
		
	}
	
	ro <- unlist(rlist)
	co <- unlist(clist)
	ret <- x[ro,co]
	attr(ret,"orders") <- list(ro,co)
	attr(ret,"ro") <- rlist
	attr(ret,"co") <- clist
	attr(ret,"nsplit") <- ncl-1
	return(ret)
	
} 
