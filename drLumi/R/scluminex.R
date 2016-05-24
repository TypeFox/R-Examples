#' Estimates different models for each analyte.
#' 
#' @description
#' Given a dilutions points and background \code{data.frame} estimates 
#' a model (in a recursive way is possible)
#' for a background method.
#'
#' @usage 
#' scluminex(plateid, standard, background, lfct, 
#'     bkg = c("ignore","subtract", "include", "constraint"), 
#'     neill.method = "finest", fmfi = "median", fec = "ec", 
#'     fanalyte = "analyte", fwell = "well", fflag = "flag", 
#'     verbose = TRUE, ...)
#'
#' @param plateid character to identify the plate
#' @param standard a \code{data.frame} with standard values information
#' @param background a \code{data.frame} with the value of the blank controls.
#' @param lfct a character vector of SelfStarting models for background method.  
#' They will be used in order if no convergence is achieved, ie: the first 
#' \code{lfct} first, if no convergence the second function, etc. Options
#' are \code{SSl5}, \code{SSl4} and \code{SSexp}.
#' @param bkg character vector specifying how the background values are 
#' treated. Options are 'ignore', 'subtract', 'include' or 'constraint'.
#' @param neill.method character specifying the grouping method for the 
#' Neill test. Default 'finest'. Other options 'c-finest', 'percentiles'
#' or the name of the grouping variable.
#' @param fmfi name of the column with MFI values
#' @param fec  name of the column with the concentration
#' @param fanalyte name of the column with the name of the analyte
#' @param fflag name of the variable with the flag to not include a 
#' record in the standard curve estimation
#' @param fwell name of the variable with the well information
#' @param verbose logical whether show the process of estimation.
#' @param ... other parameters for the model
#' 
#' @details
#' The models are fitted by the \code{nlsLM} function from the 
#' \code{minpack.lm} package. The background data can be ignore, or use to 
#' subtract the values of all MFI or be included as a point in the 
#' standard curve with a value half of the lower value of the standard points. 
#' If two or more blank controls are specified the geometric mean of the MFI 
#' is used. The names on the two datasets need to be the same and are 
#' specified by the fmfi, fec and fanalyte arguments of the function. The 
#' routine should receive the values of the MFI from the luminex fluorescence 
#' data. Analysis is performed in logarithm scale (base 10) both for the MFI 
#' and the concentrations.
#' 
#' The grouping variable for the \code{neill.method} can specified if there
#' are replicates of doses in the assay. If there are no replicates
#' one of the three 'grouping' methods can be selected.
#' 
#' @return A list with the following components model, convergence, 
#' coef, data, rsquare
#' \itemize{
#' \item{\code{model}}{, the nls model} 
#' \item{\code{convergence}}{, convergence of the model}  
#' \item{\code{coef}}{, coefficients values for the \code{nls} model}
#' \item{\code{data}}{, data of the model}
#' \item{\code{rsquare}}{, R^2 values for the performed models}  
#' }
#' 
#' @import minpack.lm  
#' @importFrom plyr ldply dlply
#' @importFrom stringr str_trim
#' @importFrom stats getInitial
#' 
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)$plate_1
#'
#' # Fit model and summary object
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#'                 lfct=c("SSl4", "SSl5"), 
#'                 bkg="ignore", 
#'                 fmfi="mfi", 
#'                 verbose=FALSE)               
#' ss <- summary(igmodels)
#' 
#' # Information
#' names(igmodels)
#' names(igmodels$FGF)
#' 
#' # Summary data
#' ss
#' as.data.frame(ss)  
#' as.data.frame(igmodels)
#' 
#' # Plot the standard curve
#' plot(igmodels,"sc")
#' 
#' 
#' @export
scluminex <- function(plateid, standard, background, lfct, 
                bkg = c("ignore","subtract","include","constraint"),
                neill.method="finest", fmfi="median", fec="ec", 
                fanalyte="analyte", fwell = "well", fflag="flag", 
                 verbose = TRUE, ...){

    if(is.null(background)){
        background <- as.data.frame(t(rep(NA,3)))
        names(background) <- c(fmfi, fec, fanalyte) 
        alertbkg <- 1
    } else {
        if(!inherits(background,"data.frame")){
            stop("'background' argument must be a data.frame")
        } else {
            if(nrow(background)==0){
                background <- as.data.frame(t(rep(NA,3)))
                names(background) <- c(fmfi, fec, fanalyte) 
                alertbkg <- 1
            } else {
                smes <- "'background' must contain mfi, ec, analyte and well."
                if(any(c(fmfi, fec, fwell, fanalyte)%nin%names(background))){
                    stop(smes)   
                } 
                alertbkg <- 0
            }       
        }  
    }
    
    if("..constraint.value"%in%c(fmfi, fec, fwell, fanalyte, neill.method)){
      stop("'..constraint.value' is a name reserved by the function" )
    }  
    
    if(!inherits(lfct,"character")) stop("'lfct' must be character")
    if(!all(lfct%in%c("SSl4","SSl5","SSexp"))){
        stop("'lfct' must be SSl4, SSl5 or SSexp")  
    } 
    if(length(unique(lfct))!=length(lfct)) stop("Specify the model only once.")

    if(!inherits(standard,"data.frame")){
        stop("'standard' argument must be a data.frame")
    } else { 
        if(nrow(standard)==0) stop("'standard' argument must have information")
        if(any(c(fmfi, fec, fwell, fanalyte)%nin%names(standard))){
            stop("'standard'must contain mfi, ec, analyte and well fields.")
        } 
    }  
    
    grouping <- NULL
    if(neill.method%in%c(fmfi, fec, fwell, fanalyte)){
      stop("'neill.method' must be different than 'fvars' args.")
    }    
    if(neill.method%in%names(standard)){
      grouping <- standard[,neill.method]
    } else{
      if(!(neill.method%in%c("finest","c-finest","percentiles"))){
        stop("neill.method argument not found")  
      }  
    } 

    la <- as.character(unique(standard[,fanalyte]))
    bkg <- match.arg(bkg)
    if(alertbkg==0)  meanbkg <- dlply(background,fanalyte,
                                function(x){geom_mean(x[,fmfi], na.rm=TRUE)})
    if( all(is.na(background[,fmfi]) ) ) alertbkg <- 1

    if(fflag%in%names(background) & nrow(background)>0){ 
        background[,fflag] <- ifelse(is.na(background[,fflag]),"",background[,fflag])
        selec_back <- background[,fflag]==""
        background <- background[selec_back,]
    }
    if(fflag%nin%names(background) & nrow(background)>0){
        background[, fflag] <- ""
    }

    if(fflag%nin%names(standard)){ 
        standard[, fflag] <- ""
    } else {
      standard[,fflag] <- ifelse(is.na(standard[,fflag]),"",standard[,fflag])
    }

    if(alertbkg==1 & bkg!="ignore"){
        stop( paste("No background recognised. No possible to '", 
                bkg, "' background only 'ignore'.", sep=""))
    } 

    if(bkg=="subtract") {
        for(s in la){
            stanalyte <- standard[which(standard[,fanalyte]==s),fmfi]
            sub.back <- stanalyte-meanbkg[[s]]
            standard[which(standard[,fanalyte]==s),fmfi] <- sub.back
        }
        if( any(standard[,fmfi]<=0, na.rm=TRUE) ) {
            standard[which(standard[,fmfi]<=0),fmfi] <- NA
        }
    }

    if(bkg=="include") {
        minbkg <- dlply(standard,fanalyte,function(x){ min(x[,fec])/2} )
        dfmfi <- ldply(meanbkg,rbind)
        names(dfmfi) <- c(fanalyte,fmfi)
        dfec <- ldply(minbkg,rbind)
        names(dfec) <- c(fanalyte,fec)
        newbkg <- merge(dfmfi,dfec,fanalyte)
        newbkg[, fwell] <- "back"
        newbkg[, fflag] <- ""
        standard <- rbind.fill(standard,newbkg)
    }

    standard[,"log10_mfi"] <- log10(standard[,fmfi] )
    standard[,"log10_concentration"] <- log10(standard[,fec])
    mod <- list()
    if(verbose==TRUE) cat(paste(plateid,":\n",sep=""))
    for(s in la) {
        if(verbose==TRUE) cat(s)
        aux_standard <- standard[which(standard[,fanalyte]==s & 
                standard[,fflag] == "" & 
                !is.na(standard[,fmfi]) & 
                !is.na(standard[,fec])) , 
                c(fanalyte,fmfi,fec,fwell,"log10_mfi","log10_concentration")]
        names(aux_standard) <- c(fanalyte,fmfi,fec,fwell,
                                "log10_mfi","log10_concentration")
    
    
    mod[[s]]$convergence <- 0
    nfunction <- 0
    niter <- length(lfct)
    
    while (mod[[s]]$convergence!=1 &  nfunction < niter ){ 
    nfunction <- nfunction + 1
    model.name <- lfct[nfunction]    
    
    if(model.name=="SSl5") n.param <- 5        
    if(model.name=="SSl4") n.param <- 4
    if(model.name=="SSexp") n.param <- 2
    
    if(bkg=="constraint"){
        cons.bkg <- log10(meanbkg[[s]])       
        n.param <- n.param - 1        
        parameters <- paste(c(letters[-1][1:n.param]), collapse=",")
        aux_standard$..constraint.value <- cons.bkg
        if(model.name=="SSl5") fct <- "SSl5cons"
        if(model.name=="SSl4") fct <- "SSl4cons"
        if(model.name=="SSexp") fct <- "SSexpcons"
        
        form <- paste0("log10_mfi ~", fct,
                       "(..constraint.value, log10_concentration,", 
                       parameters,")" )

        
        } else {
            parameters <- paste(c(letters[-1][1:n.param]), collapse=",")
            form <- paste0("log10_mfi ~", model.name,
                    "(log10_concentration," , parameters,")" )
        } 
    ## in order to be compatible with 'drc'
    st <- stats::getInitial(as.formula(form), data = aux_standard) 
    model <- suppressWarnings(try(nlsLM( form, data = aux_standard, start = st, ...), 
                 silent=TRUE))
    
        if(inherits(model,"try-error")){
            mod[[s]]$convergence <- 3
        } else {
            mod[[s]]$convergence <- get_convergence(model)
        }
    }
    mod[[s]]$data <- aux_standard
    if (mod[[s]]$convergence%in%c(2,3)){ 
        mod[[s]]$model <- model
        mod[[s]]$data$warning <- "SC_Error"
        if(mod[[s]]$convergence==2){
            mod[[s]]$data$warning <- "SC_Fit_But_No_Convergence"
        } 
        mod[[s]]$data$plateid <- plateid
        mod[[s]]$fct <- model.name
        if(alertbkg==0){
            mod[[s]]$bkg_mean <- meanbkg[[s]]
        } else {
            mod[[s]]$bkg_mean <- NA
        }  
        mod[[s]]$alertbkg <- alertbkg
        mod[[s]]$bkg_method <- bkg
    }

    if (mod[[s]]$convergence==1) {
        mod[[s]]$model <- model 
        mod[[s]]$coefficients <- summary(model)$coefficients
        mod[[s]]$rsquare <- get_rsq(model, adjusted=TRUE)
        mod[[s]]$aic <- get_aic(model)      
        mfit <- get_modelfit(model, sort(aux_standard[,"log10_concentration"], decreasing=TRUE), 
                neill.method, grouping)
        mod[[s]]$modelfit <- mfit
        mod[[s]]$neill.method <- neill.method
        mod[[s]]$data$warning <- ""
        mod[[s]]$data$predicted.log10_mfi <- predict(model)
        mod[[s]]$data$residuals <- residuals(model, type = "pearson")
        ed_mat <- suppressWarnings(get_fitted(model, bkg, meanbkg[[s]], model.name))
        mod[[s]]$data$log10.fitted.conc <- ed_mat[,1]
        mod[[s]]$data$log10.fitted.conc.se <- ed_mat[,2]
        wh <- which(is.na(mod[[s]]$data$log10.fitted.conc))
        if (length(wh)> 0) {
            wmes1 <- mod[[s]]$data[wh,"warning"]
            wmes2 <- "Not_Estimated_Concentration"
            mod[[s]]$data[wh,"warning"] <- str_trim(paste(wmes1,wmes2))
        }
        whflag <- which(standard[,fflag]!="" & standard[,fanalyte]==s)
        vars <- c(fanalyte,fmfi,fec,fwell,"log10_mfi","log10_concentration")
        mod[[s]]$flag_data <- standard[whflag, vars] 
    }
    mod[[s]]$data$plateid <- plateid
    mod[[s]]$fct <- model.name
    if(alertbkg==0){
        mod[[s]]$bkg_mean <- meanbkg[[s]]
    } else {
        mod[[s]]$bkg_mean <- NA
    }  
    mod[[s]]$alertbkg <- alertbkg
    mod[[s]]$bkg_method <- bkg
    mod[[s]]$fieldnames <- list(fmfi = fmfi, 
                        fec = fec, 
                        fanalyte = fanalyte, 
                        fwell = fwell, 
                        fflag = fflag)

    if(verbose==TRUE){
        conv.code <- mod[[s]]$convergence
        if(conv.code==1) cat(": convergence \n")
        if(conv.code==2) cat(": fit but no convergence \n")
        if(conv.code==3) cat(": error, no fit model \n")
    }  
    }
    class(mod) <- "scluminex"
    return(mod)
}
