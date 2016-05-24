sing.im <- function(data,...){
	res <- singim(data)
	return(res)
}


singim <-
function(data,...){
		n <- dim(data)[2]
		r <- dim(data)[1]	
		for(i in 1:n){
			ok <- complete.cases(data[,i])
			if(!(all(ok))){
				level <- union(data[,i],NULL)
  				id.na = which(is.na(level))
  				if(length(id.na)==1)level <- sort(level[-id.na])
				else level <- sort(level)
  				num <- length(level)
				pro <- numeric()
				for(k in 1:num)pro[k] <- length(which(data[,i]==level[k]))
				pro <- pro/sum(pro)
				miss <- which(ok==FALSE)
				for(m in miss){
					p <- runif(1)
					minus <- cumsum(pro)-p
					id <- which(minus >= 0)[1]
					data[m,i] <- level[id]
				}
			}
		}
		return(data)
}
