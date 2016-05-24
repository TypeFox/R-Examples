nearestBC_Vals <- function(framevals,WIN){
    wincoords <- data.frame(cbind(vertices(WIN)$x,vertices(WIN)$y))
    znear <- c()
    for(i in 1:nrow(wincoords)){
        tmps.x <- framevals[,1]-wincoords[i,1]
        tmps.y <- framevals[,2]-wincoords[i,2]
        hyps <- sqrt(tmps.x^2+tmps.y^2)
        hyps[is.na(framevals[,3])] <- NA
        znear <- append(znear,framevals[,3][which(hyps[!is.na(hyps)]==min(hyps[!is.na(hyps)]))])
    }
	print(znear)
    return(znear)
}
