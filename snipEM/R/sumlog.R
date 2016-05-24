
# sumlog=function(x,lower=-745,upper=709) {
  # n=length(x)
  # # log(a + b) = log(a) + log (1 + exp(log(b) - log(a)))

  # index=rep(1,n)
    # i=1
  # while(i<n & (x[i]<lower | x[i]>upper)) {
    # i=i+1}
    # if(i>n) {stop("Not Summable!")}
    # if(i<=n) {s=x[i]
             # index[i]=0}
  
  # while(sum(index)>0) {
    # idx=which(index==1)
    # l=length(idx)
    # i = 1 
 # while(i<l & ((x[idx[i]]-s)<lower | (x[idx[i]]-s)>upper)) {
    # i=i+1}
    # if(i>l) {stop("Not Summable!")}
    # if(i<=l) {s=s+log (1 + exp(x[idx[i]]-s))
             # index[idx[i]]=0}
  # }
# s}

sumlog <- function(x,lower=-745,upper=709){
	if( missing(x) ) stop("'x' missing")
	s <- tryCatch( .Call("fast_sumlog", x, lower, upper, length(x)),
		"std::range_error" = function(e){
		conditionMessage( e ) } )
	if( !is.finite(s) ) stop("Not Summable!")
	return( s )
}

