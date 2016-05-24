#' Data quality control
#'
#' This function checks data quality of anchor images by comparing the pairwise Euclidean distances of landmarks on
#' the left and right anchors for both dorsal and ventral anchors. 
#' @param x a numeric vector generated from \code{pwdist}
#' @param id.out if TRUE, allows the user to query data indices in the TMD plot
#' @param tit title for the plot
#' @param bound magnitude of upper and lower bounds on the y-axis
#' @return A matrix with three columns:
#' \item{sqasr}{the square root of average squared residuals}  
#' \item{slope}{slope of regressing M against A, multiplied by 100 to be on the same scale as average squared residuals} 
#' \item{Q}{quality score, which is given by \emph{Q=100[10^(-sqrt(sqasr^2+slope^2)/10)]} }
#' @details The quality control method implemented here is based on consideration of the difference (M) and average (A) of pairwise Euclidean distances between the left and right
#' anchors. Let the residual be defined as the deviation of M from 0. Assuming that left and right anchors are more or less symmetric, 
#' good quality data should have squared residuals that are small on average and show more or less random deviation from 0 
#' along A. The latter corresponds to slope that is near 0 in a linear regression of M against A. The result is visualized 
#' using the Tukey Mean-Difference (TMD) plot (also known as the Bland-Altman plot; see Bland & Altman (1986)), in which mean M and the 95 \% 
#' limits of agreement (within 2 standard deviations from mean M) are indicated as dashed horizontal lines. 
#' Good quality data have small values of \code{sqasr} and \code{slope}, resulting in high \code{Q} scores, and vice versa.
#' @seealso \code{\link{pwdist}}, \code{\link{boxplotSort}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Bland JM, Altman DG. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. Lancet 327:307-310.
#'
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' library(cluster)
#'
#' data(ligophorus_tpsdata)
#' data(spcolmap)
#'
#' #A low quality specimen
#' Qscore(pwdist(ligophorus_tpsdata$bantingensis[[5]],average=FALSE))
#' #A high quality specimen
#' Qscore(pwdist(ligophorus_tpsdata$bantingensis[[11]],average=FALSE))
#'
#' #Useful diagnostic plots
#' Qmat <- vector("list",length(ligophorus_tpsdata)) 
#' for(i in 1:13){
#' Qmat[[i]] <- do.call(rbind,lapply(ligophorus_tpsdata[[i]], 
#' function(k) Qscore(pwdist(k, average=FALSE))))
#' rownames(Qmat[[i]]) <- rep(names(ligophorus_tpsdata)[i],nrow(Qmat[[i]]))
#' }
#' names(Qmat) <- names(ligophorus_tpsdata)
#'
#' #Box plot for quality score by species, sorted using descending median quality score
#' Q <- lapply(Qmat, function(k) k[,3])
#' boxplotSort(Q, italic=TRUE, ylab="Quality score", df=1)
#'

Qscore <- function(x, id.out=FALSE, tit="", bound=50) {

avg <- (x[,1]+x[,2])/2
minus <- x[,1]-x[,2]

plot(avg, minus, pch=16, main=tit, xlab="A", ylab="M", ylim=c(-bound, bound))

dev <- mean(minus)
sd_dev <- sd(minus)
abline(h=c(dev, dev+2*sd_dev, dev-2*sd_dev), lty=2)

if(id.out == TRUE){
identify(avg, minus)
}

residual <- sqrt(sum(minus^2)/nrow(x))
slope <- (lm(minus~avg))$coefficients[2] * 100
delta <- sqrt(residual^2 + slope^2)
qscore <- 10^(-delta/10) * 100
result <- c(residual,slope,qscore)
names(result) <- c("sqsr","slope","Q") 
return(result)
}