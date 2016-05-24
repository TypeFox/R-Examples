#http://stackoverflow.com/questions/3789968/generate-a-list-of-primes-in-r-up-to-a-certain-number
prettyGraphsColorSelection <- function(n.colors=1,offset=NULL,starting.color=163){
	if(is.null(offset)){
		offset <- 19
	}
	##stolen
	primest <- function(n){
	    p <- 2:n
	    i <- 1
	    while (p[i] <= sqrt(n)) {
	        p <-  p[p %% p[i] != 0 | p==p[i]]
	        i <- i+1
	    }
	    return(p[length(p)])
	}	
	##stolen
	getPrime <- function(offset){
		div <- 2:floor(sqrt(offset))
		if(!any(offset %% div == 0)) {return(offset)}
		else{primest(offset)}				
	}	
	offset<-getPrime(offset)
	the.colors <- prettyGraphsColors()
	if(round(starting.color) < 0 || round(starting.color) > length(the.colors)){
		starting.color <- 163
	}	
	#the.seq <- round(seq(starting.color,length(the.colors)*n.colors,offset) %% length(the.colors))
	return(as.matrix(the.colors[(seq(1,n.colors*offset,offset) + (starting.color-1)) %% length(the.colors)]))
	
}