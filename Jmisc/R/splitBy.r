#' Split a vector by a sequence of length
#' This function will split the vector \code{x} into \code{length(x)} subvector. The length of each subvector is given by \code{by}.
#' @name splitBy
#' @aliases splitBy
#' @title Split a vector by a sequence of length  
#' @param x A vector to be splitted 
#' @param by A vector of length 
#' @return a list of subvector
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples
#' splitBy((1:10)*10,c(2,2))
#' splitBy((1:10)*10,c(2,3,4))
#' \dontrun{   
#' expect_equivalent(splitBy((1:10)*10,c(2,2)) ,  list(c(10,20),c(30,40)))
#' expect_equivalent(splitBy((1:10)*10,c(2,3,4)) , list( c(10,20), c(30,40,50) ,c(60,70,80,90)  ))
#' }

# splitBy<-function(x,by){
# 	extra<-length(x)-sum(by)
# 	if (extra<0)
# 		stop("length(x) is smaller than sum(by)")

# 	if (extra>0)
# 		x<-head(x,sum(by))

# 	index_0<-c(rep(1:length(by),times=by) )

# 	split(x,index_0)
# }
splitBy = function(x,by){
  start = c(0,head(cumsum(by),-1)) +1
  end = cumsum(by)
  out = list()
  for (i in 1:length(by)){
    if (is.matrix(x))
      out[[i]] = x[start[i]:end[i],]
    else 
      out[[i]] = x[start[i]:end[i]]
  }
  out
}

# splitBy<-function(x,by){
# 	if (length(x)<sum(by))
# 		stop("length(x) is smaller than sum(by)")
 
# 	index <- 1:by[1]
# 	out <- list(x[index])

# 	if (length(by)==1)
# 		out
# 	else
# 		c(out,splitBy(x[-index],by[-1]))
# }