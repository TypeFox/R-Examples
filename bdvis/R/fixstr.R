#'Change structure of the data frame according to the package's needs
#'
#'This function is slated to depriciate in next version. Please use function
#'\link{format_bdvis} insted.
#'
#'Modify the name of certain fields in the provided data.frame to meet the
#'requirements of the package, to make sure functions are executed properly.
#'
#'@param indf input data frame containing biodiversity data set
#'@param Latitude name of the Latitude field in original data frame
#'@param Longitude name of the Longitude field in original data frame
#'@param DateCollected name of the Date Collected field in original data frame
#'@param datefmt format string for the original date field \link[base]{strptime}
#'@param SciName name of the Scientific Name field in original data frame
#'@export
#'@family Data preparation functions
#'@examples \dontrun{
#'inat = fixstr(inat, DateCollected = "Date.collected", datefmt = "%Y-%m-%d %H:%M:%S")
#'}
fixstr <- function (indf, Latitude =NA, Longitude=NA, DateCollected=NA,
                    datefmt=NA,SciName=NA){
  .Deprecated("format_bdvis")
  chgname <- function (orig,new){
    if(is.element(orig,names(indf))){
      names(indf)[which(names(indf)==orig)]=new
      return(indf)
    } else {
      stop(paste("\n Field ",orig," not present in input data\n"))
    }
  }
  if (!is.na(Latitude)) {
    indf=chgname(Latitude,"Latitude")
  } 
  if (!is.na(Longitude)) {
    indf=chgname(Longitude,"Longitude")
  } 
  if (!is.na(DateCollected)) {
    indf=chgname(DateCollected,"Date_collected")
  } 
  if (!is.na(datefmt)) {
    indf$Date_collected=as.Date(indf$Date_collected,datefmt)
  } 
  if (!is.na(SciName)) {
    indf=chgname(SciName,"Scientific_name")
  } 
  return(indf)
}
