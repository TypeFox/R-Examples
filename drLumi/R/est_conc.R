#' Estimates concentration given an scluminex object
#' 
#' @description
#' Given a \code{scluminex} object with standard curve information 
#' and a \code{data.frame} with response values add to the original 
#' dataset the concentration data.
#' 
#' @usage
#' est_conc(x, df, fanalyte = "analyte", fmfi = "median", dilution = 1,
#'     one.curve = FALSE, level = 0.95)
#' 
#' @param x a \code{scluminex} object.
#' @param df input \code{data.frame} with the analyte and median 
#' fluorescence intensity variables.
#' @param fanalyte name of the field with the analyte information. 
#' Default 'analyte'.
#' @param fmfi name of the field with the mfi (response) information. 
#' Default 'median'.
#' @param dilution numeric value of the dilution that must be used for 
#' the estimation of the concentration.
#' @param one.curve logical according if only one curve must be used for 
#' estimation.
#' @param level numeric value, confidence level, required. Default 0.95.
#'  
#' @details
#' Given a \code{\link{scluminex}} object and a \code{data.frame} with 
#' analyte and MFI information adds the concentration information to the 
#' dataset (concentration, standard error of the concentration and a 
#' warning variable). The MFI data will be transformed into log10(MFI). 
#' The method utilized is the Delta method of \code{\link{invest}} function.
#' 
#' Merging variables are defined in the \code{fanalyte} and \code{fmi} 
#' arguments of the function.
#' 
#' If only one standard curve is fitted for several analytes \code{one.curve} 
#' argument must be specified to \code{TRUE} and \code{scluminex} must have
#' only one analyte information. The same \code{scluminex} information will
#' be used for all analytes of the \code{df} \code{data.frame}.
#'  
#' If one standard curve is fitted by each analyte \code{one.curve} must 
#' be \code{FALSE}, so the function will merge each model of the \code{scluminex}
#' object with the corresponding analyte of the \code{df} argument.
#' 
#' 
#' @return
#' Input \code{data.frame} with the following merged variables:
#' \itemize{
#' \item{\code{log10.fitted.conc}}{, log10 fitted concentration} 
#' \item{\code{log10.fitted.conc.se}}{, log10 standard error of the log10 
#' fitted concentration}  
#' \item{\code{dilution}}{, dilution to be applied to the samples}
#' \item{\code{dil.fitted.conc}}{, diluted fitted concentration 
#' in original scale}
#' \item{\code{dil.lb.conc}}{, diluted fitted lower bound concentration 
#' in original scale}  
#' \item{\code{dil.ub.conc}}{, diluted fitted upper bound concentration 
#' in original scale}  
#' \item{\code{warning}}{, warning message (if necessary)}  
#' }
#'  
#'  @examples
#' # Load data and fit models  
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' sdf <- data_selection(dat, ecdata)$plate_1
#' 
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background,
#'             lfct="SSl4", bkg="ignore", fmfi="mfi", verbose=FALSE)
#'             
#' # Data to estimate concentration
#' concdf <- sdf$positive
#'       
#' # Dilution factor 1
#' est_conc(igmodels, concdf, fmfi="mfi", dilution = 1) 
#' 
#' # Dilution factor 2
#' est_conc(igmodels, concdf, fmfi="mfi", dilution = 2)  
#' 
#' 
#'  @importFrom stringr str_trim
#'  @importFrom plyr ldply dlply rbind.fill ddply
#'  @export
est_conc <- function(x, df, fanalyte="analyte", fmfi="median", 
                dilution=1, one.curve = FALSE, level = 0.95){

    stopifnot(inherits(x,"scluminex"))

    if(length(x)>1 & one.curve==TRUE){
        stop("More than 1 model estimated but one.curve set to TRUE")
    }
    if(!inherits(one.curve, "logical")) stop("'one.curve'  must be logical")
    if(!inherits(dilution, "numeric")) stop("'dilution'  must be numeric")
    if(!inherits(df, "data.frame")) stop("'df' must be a data.frame")
    if(level<0 | level>1) stop("'level'  must be a value between 0 and 1")  
    if(fanalyte%nin%names(df)) stop("'fanalyte' not found in 'df'")
    if(fmfi%nin%names(df)) stop("'fmfi' not found in 'df'")

    ans <- list()
    numrows <- nrow(df)
    la <- unique(as.character(df[,fanalyte]))

    for (s in la) {
        if(one.curve) names(x) <- s
        ans[[s]] <- df[which(df[,fanalyte]==s),]
        if (x[[s]][["bkg_method"]]== "subtract"){
            ans[[s]][,fmfi] <- ans[[s]][,fmfi]-x[[s]][["bkg_mean"]]
        }
        ans[[s]][,"warning"] <- ""
        if(x[[s]]$convergence==1){
            if (s %in% names(x)) {
                model <- x[[s]]$model
                log10mfi <- log10(ans[[s]][,fmfi])
                fct <- x[[s]]$fct 
                bkg_mean <- NULL
                if(x[[s]][["bkg_method"]]=="constraint"){
                    bkg_mean <- x[[s]][["bkg_mean"]]
                }
            ecdf <- invest(x, s, log10mfi, "delta", level)

            dil.fitted.conc <- 10^(ecdf[,2] + log10(dilution))
            dil.lb.conc <- 10^(ecdf[,3] + log10(dilution))
            dil.ub.conc <- 10^(ecdf[,4] + log10(dilution))
            newinfo <- cbind(ecdf[,2],ecdf[,5], dilution, dil.fitted.conc, 
                        dil.lb.conc, dil.ub.conc)
            newinfo <- as.data.frame(newinfo)
            ori.vars <- c("log10.fitted.conc", "log10.fitted.conc.se")
            dil.vars <- c("dil.fitted.conc","dil.lb.conc","dil.ub.conc")
            names(newinfo) <- c(ori.vars,"dilution",dil.vars)

            if(any(names(newinfo)%in%names(ans[[s]]))){
                stop("Some of the new added variables are also in the 'df'")    
            }
            ans[[s]] <- cbind(ans[[s]],newinfo)

            ## Warnings for concentrations not estimated  or beyond limits
            ans[[s]][,"warning"] <- ""
            wh <- which(is.na(ans[[s]][,"log10.fitted.conc"]))
            war.mes <- str_trim(paste(ans[[s]][wh, "warning"],"NOT_ESTIMATED"))
            ans[[s]][wh, "warning"] <- war.mes
            } else {
                ans[[s]][, "warning"] <- "ANALYTE_NOT_FOUND_IN_SCLUMINEX"
                ans[[s]][,"log10.fitted.conc"]<-NA
                ans[[s]][,"log10.fitted.conc.se"]<-NA
            }
        } else {
            ans[[s]][,"warning"]<-"SC_NON_CONVERGENCE"
            ans[[s]][,"log10.fitted.conc"]<-NA
            ans[[s]][,"log10.fitted.conc.se"]<-NA
        }
    }
    ans <- ldply(ans,rbind)
    ans$.id <- NULL
    stopifnot(nrow(ans)==numrows)
    return(ans)
}





