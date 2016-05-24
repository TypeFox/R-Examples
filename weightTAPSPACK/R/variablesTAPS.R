#' Lists names of variables by wave
#'
#' Lists the names of the variables in the TAPS dataset for the wave(s) specified
#'
#' @param month A character vector specifying the first three letters of the month(s) of the wave(s) of interest
#' @param year A numeric vector specifying the year(s) of the wave(s) desired. This vector MUST BE the same length of the month vector (one year per month specified)
#'
#' @return A list with the names of the variables per wave specified
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu}  Michelle Torres: \email{smtorres@@wustl}
#' @rdname variablesTAPS
#' @aliases variablesTAPS,ANY-method
#' @examples
#' 
#' variablesTAPS(month=c("Feb","Mar"), year=c(2012,2012))
#' @seealso \code{\link{weightTAPS}} \code{\link{weightTAPSPACK}} \code{\link{subsetTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @export
setGeneric(name="variablesTAPS",
           def=function(month, year)
           {standardGeneric("variablesTAPS")}
)

setMethod(f="variablesTAPS",
          definition=function(month, year){
  month_year <- data.frame(Year=rep(c(2012,2013),each=12),
                           Month=rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),2),
                           Wave=2:25)
  num.waves <- paste("S", 
                     apply(cbind(month,year),1,function(x) month_year$Wave[month_year$Month==x[1] & month_year$Year==x[2]]), 
                     sep="")
  scan.wave <- function(word){
    l.word <- length(word)
    pos.S <- tail(grep("S", word),1)
    term <-ifelse(length(pos.S)==0, "NW", paste(word[pos.S:l.word], collapse=""))
    return(term)
  }
  wave.end <- laply(strsplit(colnames(TAPScum), split=""), function(x) scan.wave(x))
  selected.columns <- sapply(num.waves, function(z) grep(pattern=z, x=wave.end))
  if(length(num.waves)>1){
    selected.variables <- llply(selected.columns, function(x) colnames(TAPScum)[x])
  }
  else{selected.variables <- list(colnames(TAPScum)[as.numeric(selected.columns)])}
  names(selected.variables) <- paste(month, year, sep=" ")
  return(selected.variables)
}
)
