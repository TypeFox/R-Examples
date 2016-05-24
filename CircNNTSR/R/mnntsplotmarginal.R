mnntsplotmarginal<-function(cestimates,M,component, ...){

    nntsplotint <- function(theta) {
        res <- mnntsmarginal(cestimates,M,component,theta)
	return(res)
    }
    return(curve(nntsplotint, 0, 2 * pi, ylab="Marginal density function",xlab = paste("Component:",component), ...))
}