#' Makes a list used in \code{pointkriging}
#' 
#' Makes a list of four components sequences based on input number for use in
#' \code{pointkriging}.
#' 
#' 
#' @param nm A number
#' @return Returns a list with components: \item{rrt}{sequence 1 of length
#' \code{4*nm^2}} \item{crt}{sequence 2, same length as \code{rrt}}
#' \item{rrt}{sequence 3, same length as \code{rrt}} \item{crt}{sequence 4 of
#' length \code{nm + 1}} for use in \code{pointkriging}.
#' @note Needs further elaboration.
#' @seealso \code{\link{pointkriging}}
#' @keywords manip
#' @export bua
bua <-
function(nm = 10)
{
	rrt <- c(0, 0, 1, 1)
	crt <- c(0, 1, 1, 0)
	for(i in 2:nm) {
		stdrrt <- c(matrix(0, 8 * i - 4, 1))
		stdcrt <- stdrrt
		n <- i * 8 - 4
		stdrrt[1:(2 * i)] <- 1 - i
		stdrrt[(2 * i + 1):(4 * i - 1)] <- c((2 - i):i)
		stdrrt[(4 * i - 1):(6 * i - 2)] <- i
		stdrrt[(6 * i - 2):(8 * i - 4)] <- c(i:(2 - i))
		stdcrt[1:(2 * i)] <- (c(1 - i):i)
		stdcrt[(2 * i + 1):(4 * i - 1)] <- i
		stdcrt[(4 * i - 1):(6 * i - 2)] <- c(i:(1 - i))
		stdcrt[(6 * i - 2):(8 * i - 4)] <- 1 - i
		crt <- c(crt, stdcrt)
		rrt <- c(rrt, stdrrt)
	}
	i1 <- 4
	for(i in 2:nm) {
		i1[i] <- i1[i - 1] + 8 * i - 4
	}
	i1 <- c(0, i1)
	# 	Part that comes instead of the loop that
	#	is too slow.  
	ind <- c(3, 0, 2, 0, 4, 0, 1)
	r1 <- rrt - 0.1
	c1 <- crt - 0.1
	rr <- sign(r1) + sign(c1) * 2 + 4
	dir <- ind[rr]
	#	dir<- c(1:length(rrt))
	#	for(i in 1:length(rrt)){
	#		if(rrt[i]>0 && crt[i]>0)dir[i]<-1
	#		if(rrt[i]<= 0 && crt[i]>0)dir[i]<-4
	#		if(rrt[i]<= 0 && crt[i]<=0 )dir[i]<-3
	#		if(rrt[i]>0 && crt[i]<=0)dir[i]<-2
	#	}
	return(list(rrt = rrt, crt = crt, dir = dir, i1 = i1))
}

