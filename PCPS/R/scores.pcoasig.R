#' @rdname pcoa.sig
#' @encoding UTF-8
#' @method scores pcoasig
#' @export
scores.pcoasig<-function (x, choices=c(1,2), ...) 
{
	if (length(choices) != 2) {
        stop("\n Choices must have length equal to two \n")
    }
    max1 <- max(x$PCoA$vectors[, choices[1]])
    max2 <- max(x$PCoA$vectors[, choices[2]])
    min1 <- min(x$PCoA$vectors[, choices[1]])
    min2 <- min(x$PCoA$vectors[, choices[2]])
    scores1<-ifelse(x$correlations[,choices[1]]>0,x$correlations[,choices[1]]*max1,x$correlations[,choices[1]]*abs(min1))
	scores2<-ifelse(x$correlations[,choices[2]]>0,x$correlations[,choices[2]]*max2,x$correlations[,choices[2]]*abs(min2))
	rscores<-cbind(scores1,scores2)
	colnames(rscores)=colnames(x$vectors[, choices])
    Res <- list(scores.sites = x$PCoA$vectors[, choices], 
        scores.species = rscores)
    return(Res)
}