rebuildMCAtable <- function(DATA){
	##private function
	mat.rep <- function(fill.matrix,which.column,times,items){
		na.locs <- which(is.na(fill.matrix[,which.column]))
		rep.val <- rep(items,times)
		fill.matrix[na.locs[1:length(rep.val)],which.column] <- rep.val
		return(fill.matrix)
	}
	
	na.locs <- which(DATA < 1 & DATA > 0)
	colSums.data <- colSums(DATA)
	orig.cols <- sum(cumsum(colSums.data) %% nrow(DATA)==0)	
	if(length(na.locs)!=0){

		DATA.copy <- DATA
		DATA.copy[na.locs] <- 0
			colSums.data.copy <- colSums(DATA.copy)	
		fill.mat <- matrix(data=NA,nrow(DATA),orig.cols)
		
		actuals <- c(colSums.data-(colSums.data-colSums.data.copy))
		end.locs <- which(cumsum(colSums.data) %% nrow(DATA)==0)
		begin.locs <- c(0,end.locs[1:(length(end.locs)-1)])
		num.vars <- (end.locs-begin.locs)
	
		col.indices <- rep(1:orig.cols,num.vars)
		items <- 1:length(colSums(DATA))
		
		#fill.mat.copy <- fill.mat
		##I NEED TO BE CRAFTY HERE. I can probably use one of the apply functions...
		for(i in 1:length(col.indices)){
			fill.mat <- mat.rep(fill.mat,col.indices[i],actuals[i],items[i])
		}
		#DATA <- apply(fill.mat,2,sample)
		DATA <- fill.mat
	}else{
		DATA <- matrix(unlist(mapply(rep,times=colSums.data,x=as.matrix(1:length(colSums.data)))),nrow(DATA),orig.cols)
	}
	return(DATA)
}