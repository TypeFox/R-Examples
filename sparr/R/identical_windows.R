identical_windows <- function(w1,w2){
	w1x <- vertices(w1)$x
	w1y <- vertices(w1)$y
	w2x <- vertices(w2)$x
	w2y <- vertices(w2)$y
	
	if(!all(c(length(w1x)==length(w2x),length(w1y)==length(w2y)))) return(FALSE)
	
	if(any(c(sum(w1x!=w2x),sum(w1y!=w2y))>0)) return(FALSE)
	else return(TRUE)
}
	