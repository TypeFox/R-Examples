#' Median-ordered box plot
#'
#' This is a modification of the standard box plot to produce box plots ordered on the basis of descending median, 
#' for more effective graphical representation of data.
#' @param x a list, with species as the names of the list
#' @param italic if TRUE, the species names are italicized
#' @param srt the angle of the x-axis labels relative to the x-axis
#' @param df displacement factor for positioning x-axis labels
#' @param varwidth if TRUE, the width of the box plot for a species is proportional to the square root of the species's sample size 
#' @param col a character vector specifying the box colors for two clades
#' @param ylab y-axis title
#' @param main title for the plot
#' @param clade a character vector specifying the species names for the first clade; currently supports only two clades
#' @details The box plots produced using this function are helpful for documenting sample quality score distributions 
#' of each species. They are also useful for checking the distribution of pairwise Euclidean distance between two landmarks.
#' @seealso \code{\link{Qscore}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' data(pwed_pd)
#' pwed_pd <- matrix2list.2(pwed_pd)
#'
#' cladeI <- c("liewi","fenestrum","grandis","johorensis","kedahensis","kederai")
#' #We just want to look at distance between LM1 and LM3 in for dorsal anchor
#' boxplotSort(lapply(pwed_pd, function(k) k[,which(colnames(k)=="D1_3")]), italic=TRUE,
#' col=c("dodgerblue","violetred"), clade=cladeI, 
#' ylab=expression(paste("Length ", "(", italic(mu),"m", ")")))
#'
#' #Separation of two lineage seems possible at 15 micrometers
#' abline(h=15)
#'

boxplotSort <- function(x, italic=FALSE, srt=45, df=0.5, varwidth=TRUE, col=NULL, ylab="", main="", clade=NULL){
med <- unlist(lapply(x, median))
ord <- order(med, decreasing=TRUE)
x_ord <- vector("list",length(x))
for(i in 1:length(x_ord)){
x_ord[[i]] <- x[[ord[i]]]
}
names(x_ord) <- names(x)[ord]
boxcol <- rep("white", length(x))

if(is.null(clade) == FALSE){
cladeid <- which(names(x_ord) %in% clade) 

boxcol[cladeid] <- col[1]
boxcol[-cladeid] <- col[2]
}

boxplot(x_ord, xaxt="n",xlab="Species",ylab=ylab,varwidth=varwidth, col=boxcol, main=main, medlwd=1)

    splabel <- names(x_ord)
    
    f <- splabel
    if(italic==TRUE){
    f <- vector("expression", length(splabel))
    for (s in 1:length(splabel)) {
        f[[s]] <- substitute(italic(nn), list(nn = splabel[s]))
    }
    }

axis(1, 1:length(x), labels=FALSE)
text(1:length(x), par("usr")[3] - df, labels=f, srt=srt, pos=1, xpd=TRUE, cex=0.6)
}

