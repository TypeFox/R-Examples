#'  Extract outliers from a scluminex object
#'  
#' @description
#' This function extracts outliers based on a residuals cutoff from 
#' \code{scluminex} object.
#'  
#' @param x a \code{scluminex} object
#' @param out.limit cutoff point for outliers
#' 
#' @details The residuals are the standardized ones 
#' for \code{\link{nls}} models.
#' 
#' @return A \code{data.frame} with samples whose standardized 
#' residual is greater than \code{out.limit}.
#' 
#' @examples
#' # Load data and estimate models
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' # Select data and fit models
#' sdf <- data_selection(dat, ecdata)[[1]]
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#'                       lfct=c("SSl4", "SSl5"), 
#'                       bkg="ignore", 
#'                       fmfi="mfi", 
#'                       verbose=FALSE)
#' 
#' # Extract outliers > abs(2.5)
#' get_outliers(igmodels, out.limit=2.5)
#'
#' # Extract outliers > abs(1.5)
#' get_outliers(igmodels, out.limit=1.5)
#' 
#' # Extract outliers > abs(1)
#' get_outliers(igmodels, out.limit=1)
#'  
#'    
#' @export 
get_outliers <- function(x, out.limit=2.5) {  
    if(!inherits(x,"scluminex")) stop("'x' must be a scluminex object")
    if(!inherits(out.limit,"numeric")){
        stop("'outlier' must be a numeric object")
    }
    fwell <- x[[1]]$fieldnames$fwell
    fanalyte <- x[[1]]$fieldnames$fanalyte
    ans_data <- as.data.frame(x)
    cwells <- ans_data[,fwell]
    canalytes <- ans_data[,fanalyte]
    cplateid <- ans_data[,"plateid"]
    ans_data$batch_well_analyte <- paste(cplateid,cwells,canalytes,sep="*")
    ans_data$flag <- "OUTLIER"
    ans <- subset(ans_data,abs(residuals)>out.limit)
    vars <- c("batch_well_analyte","plateid",fwell,fanalyte,"flag","residuals")
    ans <- ans[,vars]
    nn <- c("batch_well_analyte","batch","well","analyte","flag","observations")
    names(ans) <- nn
    return(ans)
}
