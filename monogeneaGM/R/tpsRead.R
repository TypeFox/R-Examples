#' Reading tps files into R
#'
#' This function converts all tps files in a folder to a list of matrices. The tps files are created 
#' using the TPSDIG2 software (Rohlf, 2013, 2015).   
#' @param x a constant that indexes the folder of interest
#' @return a list of matrices containing the xy-coordinates of landmarks
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Rohlf FJ. (2013). Morphometrics at SUNY Stony Brook. Available at http://life.bio.sunysb.edu/morph/soft-dataacq.html.
#'
#' Rohlf FJ. (2015). The tps series of software. Hystrix. 26, doi: 10.4404/hystrix-26.1-11264.
#'
#' @examples
#' ##not run
#' ##download the tps data folders from Khang et al. (2015) and unzip them into a folder called ligotps
#'
#' #setwd("ligotps/")
#' #dir()
#' ##read tps files from the /bantingensis folder and turn them into objects of a list
#' #xydata <- tpsRead(1)
#'

tpsRead <- function(x){

	tpsfiles <- list.files(dir()[x]) 
	L <- vector("list",length(tpsfiles))
	nr <- numeric(length(tpsfiles))

	for(i in 1:length(tpsfiles)){
		rip <- scan(file=paste(c(dir()[x],"/",tpsfiles[i]), collapse=""), nlines=1, what="character")
		nr[i] <- as.numeric(unlist(strsplit(rip,"="))[2])
		L[[i]] <- scan(file=paste(c(dir()[x],"/",tpsfiles[i]), collapse=""), skip=1, nlines=nr[i])
		L[[i]] <- matrix(L[[i]],nr[i],2, byrow=TRUE)
	}

	names(L) <-  tpsfiles
	if(var(nr) > 0) {
	error <- which(nr != median(nr)) 
	warning(paste("These files:",tpsfiles[error],"have number of rows that differ from the majority of other files"))
	}
	return(L)

}