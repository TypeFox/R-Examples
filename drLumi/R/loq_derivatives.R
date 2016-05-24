#' Limits of quantifications estimation using derivatives
#' 
#' @description
#' Estimates the limits of quantification based on the second 
#' order derivative of the functions.
#' 
#' @param x a \code{scluminex} class object
#' @param subset.list list of analytes to estimate. 
#' Default \code{NULL} (all analytes of the \code{scluminex} object)
#' @param ... further arguments to be passed to \code{uniroot.all} function.
#' 
#' @details The limits of quantification are based on the 
#' maximum and minimum points of the second order derivative of the function.
#' The solution is based on the specific expressions of the derivatives. 
#' The \code{\link[rootSolve]{uniroot.all}} function is used
#' in order to find the global maximum and minimum.
#' 
#' @return Object of class \code{loq}.
#' 
#' @references 
#' Ritz C and Spiess AN (2008).
#' qpcR: an R package for sigmoidal model selection in q
#' uantitative real-time polymerase chain reaction analysis.
#' \emph{Bioinformatics} \bold{24}, 1549-51.
#' 
#' @importFrom rootSolve uniroot.all
#' 
#' @examples 
#' # Load data and estimate models  
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)[[1]]
#' 
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#'                 lfct=c("SSl4", "SSl5"), 
#'                 bkg="ignore", 
#'                 fmfi="mfi", 
#'                 verbose=FALSE)
#' 
#' loq_derivatives(igmodels)
#' 
#' @export
loq_derivatives <- function(x, subset.list=NULL, ...){
    if(!inherits(x,"scluminex")) stop("'x' must be 'scluminex' class")
    if(is.null(subset.list)){
        lanalytes <- sort(names(x))
    } else {
        if(any(subset.list %in% names(x)==FALSE)){
        stop(" 'subset.list' vector not found in 'x' ") 
    }   
        lanalytes <- subset.list
    }
    ans <- list()
    for(i in lanalytes){
        msg <- NULL
        if(x[[i]]$convergence==1){
            model <- x[[i]]$model    
            bkg.method <- x[[i]]$bkg_method
            bkg.mean <- x[[i]]$bkg_mean
            fct <- x[[i]]$fct      

            lhs <- as.character(model$m$formula()[[2]])
            parameters <- names(model$m$getPars())
            allobj <- ls(model$m$getEnv())
            rhs <- allobj[-match(c(parameters,lhs),allobj)]
            ndf <- data.frame(get(lhs, model$m$getEnv()), 
                                get(rhs, model$m$getEnv()))
            names(ndf) <- c(lhs, rhs)

            if(fct=="SSl4") FUN <- d3SSl4
            if(fct=="SSl5") FUN <- d3SSl5
            if(fct=="SSexp") FUN <- d3SSexp
  
            ff <- function(x){
                ans <- do.call(FUN, as.list(c(x,arguments)))
                ans  
            }  

            if(bkg.method!="constraint"){
                arguments <- as.numeric(coef(model))
            } else{          
                if(fct=="SSl4"){
                    arguments <- as.numeric(c(coef(model)[1], log10(bkg.mean), 
                                coef(model)[2], coef(model)[3]))
                }  
                if(fct=="SSl5"){
                    arguments <- as.numeric(c(coef(model)[1], log10(bkg.mean), 
                                coef(model)[2], coef(model)[3], 
                                coef(model)[4]))
                }  
                if(fct=="SSexp"){
                    arguments <- as.numeric(c(coef(model)[1], log10(bkg.mean)))
                }  
            }    

            ff <- Vectorize(ff)
            roots <- try(uniroot.all(ff, c(min(ndf[,rhs]),max(ndf[,rhs])), ...))
            if(inherits(roots,"try-error")){
                roots <- c(NA, NA)
                msg <- " (no estimated, error in 'uniroot.all' function)"
            } 
            if(length(roots)>2L){
                roots <- c(NA,NA)
                msg <- " (no estimated, >2 roots in third derivative)"
            } 
            if(length(roots)<2L){
                roots <- c(NA,NA)
                msg <- " (no estimated, <2 roots in third derivative)"
            }

            lloq <- min(roots)
            uloq <- max(roots)

            ly <- conf_bands(x, i, lloq)[,1]
            uy <- conf_bands(x, i, uloq)[,1]

            origdata <- ndf[,c(lhs,rhs)]

            ans[[i]] <- list(lloq=lloq, uloq=uloq, 
                        method = paste0("derivative", msg),ly=ly, uy=uy,  
                        data = origdata)
    } else {
        mes <- "(no estimated, no convergence model)"
        ans[[i]] <- list(lloq=NA, 
                    uloq=NA, 
                    method= paste0("derivatives ",mes),                     
                    uy=NA, ly=NA, ul=NA, ll=NA,  data=NA)
    }
        
    }
    class(ans) <- c("loqderivatives", "loq")  
    return(ans)
}


