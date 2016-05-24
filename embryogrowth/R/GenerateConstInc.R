#' GenerateConstInc generates a data.frame with constant incubation temperature and incubation duration
#' @title Generate a data.frame with constant incubation temperature and incubation duration
#' @author Marc Girondot
#' @return A date.frame that can be used with FormatNests()
#' @param durations A vector with incubation durations
#' @param temperatures A vector with incubation temperatures
#' @param names A vector of column names
#' @description Generate a data.frame from constant incubation temperature and incubation duration
#' @examples
#' \dontrun{
#' temp_cst <- GenerateConstInc(durations=c(150000, 100100, 100000), 
#'	temperatures=c(28, 30.5, 30.6), 
#' 	names=c("T28", "T30.5", "T30.6"))
#' }
#' @export

GenerateConstInc <- function(durations=stop("At least one incubation length must be provided"), 
                          temperatures=stop("At least one incubation temperature must be provided"), 
                          names=NULL) {
  
  if (is.null(names)) names <- paste0("X", as.character(temperatures))
  if (length(durations)!=length(temperatures)) {
    print("Same number of temperatures and incubation durations must be provided")
    return()
  }
  temp.constant <- data.frame(Time=0, as.list(temperatures))
  colnames(temp.constant) <- c("Time", names)
  temp.constant <- rbind(temp.constant, 
                            temp.constant[rep(1, length(temperatures)),])
  temp.constant[2:(length(temperatures)+1), 1] <- durations[order(durations)]
  for(i in 2:dim(temp.constant)[1]) {
    if (temp.constant[i-1, "Time"]==temp.constant[i, "Time"])
    temp.constant[i-1, 1] <- NA
  }
  temp.constant <- temp.constant[!is.na(temp.constant$Time), ]
  
  for(i in 2:(length(temperatures)+1)) {
    if (which(temp.constant$Time==durations[i-1])!=dim(temp.constant)[1])
      temp.constant[(which(temp.constant$Time==durations[i-1])+1):dim(temp.constant)[1], i] <- NA
  }
  
  return(temp.constant)
}
