#' @export 
precintcon.ci.analysis <- function(..., interval = 1, args = NA) {

  l <- list(...)
	
  if (length(l) < 1){
    stop("precintcon.ci.analysis function called without input data.");
  }
		
	set <- NULL
	
	pars <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)
	
	for (i in 1:length(l)) {
		
		data <- NULL
		
		object <- l[[i]]
		
		if (is.element("precintcon.daily", class(object)))
			object <- precintcon.fd(precintcon.classification(object, interval))
		
		if (is.element("precintcon.fd", class(object))) {
			
			a  <- exp(precintcon.ln.a(object$p.sum.n, object$p.sum.P))
			b  <- precintcon.b(object$p.sum.n, object$p.sum.P)
			r2 <- precintcon.r.squared(object, a, b * 100)
			A  <- (a/b) * (exp(1)^(b*100)) * (100 - (1/b))
			S  <- 5000 - A
			ci <- 2 * S / 10000
			
			data <- data.frame(a=a, b=b, r2=r2, A=A, S=S, ci=ci)
			
			class(data) <- c("data.frame", "precintcon.ci")
	
			if (length(l) == 1)
				return(cbind(dataset=paste(pars[i], sep=""),data))
			else
				set <- rbind(set, cbind(data.frame(dataset=paste(pars[i], sep="")), data))
			
		} else 
			stop("object should be of type \"precintcon.fd\"")
	}
		
	class(set) <- c("data.frame", "precintcon.ci")
			
	return(set)
}
