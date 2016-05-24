#' @title Point Biserial Correlations
#' @export ptbis
#' @description Calculation of the point biserial correlations for dicho- or polytomous item categories with total scale (person parameter).
#'  
#' @details no details in the moment.
#' @param y either an object of class \code{"pers"}, or an numeric vector as an result of any scaling approach (WLE, MLE, Rawscore, etc. ) relating to the Items (columns) in \code{daten}.    
#' @param daten if argument y is not an object of class \code{"pers"}, a \code{"data.frame"}, potentially with missing values, comprising dicho- or polytomous items (columns). 

#' @return An object of class \code{c("data.frame", "ptbis")} containing item statistics.
#' @exportClass ptbis
#' @examples ######################
#' ########
#' data(sim200x3) # loading reponse data
#' y <- rowSums(sim200x3)
#' ptbis(y=y, daten=sim200x3)
#' #### 
#' result <- pers(pair(sim200x3))
#' ptbis(y= result)

ptbis <- function(y, daten=NULL){# daten = data.frame mit polytomen variablen; y =skala
  
  if(any(class(y)=="pers")){
    daten <- as.data.frame(y$pair$resp)
    abil <- y$pers$WLE
  }
  if(is.vector(y)==TRUE){
    daten <- daten
    abil <- y
  }
  if(class(daten)!="data.frame"){stop("daten must be a data.frame")}
  stopifnot(dim(daten)[1]==length(abil))
  datenl <- as.list(daten)
  erglist <- lapply(datenl, polyptb, abil)
  ##### 
  maxLen <- max(sapply(erglist, length))
  # create a new list with elements padded out with NAs
  newM <- lapply(erglist, function(.ele){c(.ele, rep(NA, maxLen))[1:maxLen]})
  ergmat <- do.call(rbind, newM)
  
  colnames(ergmat) <- c(sapply(0:((dim(ergmat)[2]/2)-1),function(x){paste(x,c("rptb","n"),sep=".")}))
 
  return(ergmat)
}
