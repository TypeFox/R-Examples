#' Limits of quantifications estimation using coefficient of variation
#' 
#' Estimates the limits of quantification based on an approximation 
#' of the coefficient 
#' of variation. 
#' 
#' @details For each value of the response, the estimated concentration 
#' value and the approximated standard deviation is estimated. 
#' The function \code{\link{invest}} with the delta method approach is used. 
#' The coefficient of variation of the log10 concentration is calculated as the
#' \deqn{\sqrt{e^{ (SE \times  log(10))^2} - 1 }}
#' 
#' @param x a \code{scluminex} object
#' @param subset.list list of analytes to estimate. 
#' Default \code{NULL} (all analytes of the \code{scluminex} object).
#' @param max.cv is the target coeficient of variation by default 0.2 
#' @param n.cuts is the number of cuts to search the coefficent of variation. 
#' Default 100.
#' 
#' @return Object of class \code{loq}.
#' 
#' @references 
#' Gottschalk PG, and Dunn JR. (2005). 
#' Determining the error of dose estimates and minimum and maximum 
#' acceptable concentrations from assays with nonlinear dose-response curves. 
#' \emph{Comput Methods Programs Biomed} \bold{80}, 204-215.
#' 
#' Defawe OD, Fong Y, Vasilyeva E, Pickett M, Carter DK, Gabriel E, 
#' Rerks-Ngarm S, Nityaphan S, Frahm N, McElrath MJ and De Rosa SC.(2012).
#' Optimization and qualification of a multiplex bead array to assess cytokine 
#' and chemokine production by vaccine-specific cells.
#' \emph{J Immunol Methods} \bold{382}, 117-128.
#'  
#' 
#' @examples 
#' # Load data and estimate models  
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)$plate_1
#' 
#' igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
#' lfct=c("SSl4", "SSl5"), bkg="ignore", fmfi="mfi", verbose=FALSE)
#' 
#' loq_cv(igmodels, max.cv=0.25, n.cuts=100)
#' 
#' @export
loq_cv <- function(x, subset.list=NULL, max.cv=0.2, n.cuts = 100 ){

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
        if(x[[i]]$convergence==1){
            model <- x[[i]]$model
            bkg_method <- x[[i]]$bkg_method
            fct <- x[[i]]$fct
            bkg_mean <- x[[i]]$bkg_mean
            lhs <- as.character(model$m$formula()[[2]])
            parameters <- names(model$m$getPars())
            allobj <- ls(model$m$getEnv())
            rhs <- allobj[-match(c(parameters,lhs),allobj)]
            ndf <- data.frame(get(lhs, model$m$getEnv()), 
                        get(rhs, model$m$getEnv()))
            names(ndf) <- c(lhs, rhs)
        
            yrange <- seq(min(ndf[,lhs]),max(ndf[,lhs]),length.out = n.cuts)
            inv <- suppressWarnings(invest(x, i, yrange, "delta"))

            inv$cv <-  sqrt( exp( (inv[,5]*log(10))^2) - 1 )
            grid.data <- inv
            grid <- grid.data[!is.na(grid.data$cv),]
        
            if(min(grid$cv)<=max.cv){
                grid$inv.fit <- grid[,2]  
                grid <- grid[order(grid$inv.fit),]
                xfitmin <- grid$inv.fit[which(grid$cv==min(grid$cv))]
                loqcrit <- max(grid$cv[grid$cv<=max.cv & grid$inv.fit<=xfitmin])
                uloqcrit <- max(grid$cv[grid$cv<=max.cv & grid$inv.fit>=xfitmin])
                obslloq <- which(grid$cv==loqcrit)
                obsuloq <- which(grid$cv==uloqcrit)
    
                lloq <- max(grid$inv.fit[obslloq], min(ndf[,rhs]), na.rm=TRUE)
                uloq <- min(grid$inv.fit[obsuloq], max(ndf[,rhs]), na.rm=TRUE)
                ly <- conf_bands(x,i, lloq)[,1]
                uy <- conf_bands(x,i, uloq)[,1]
                grid$inv.fit <- NULL
            
                ans[[i]] <- list(lloq = lloq, uloq=uloq, 
                            method=paste0("coefficient variation <=", max.cv), 
                            ly=ly, uy=uy, cv=grid, data=ndf, grid.data = grid.data)
                
            } else {
              mes <- "(no estimated, minimum CV > 'max.cv')"
              ans[[i]] <- list(lloq=NA, uloq=NA, 
                               method=paste0("coefficient variation <=",
                                             max.cv, mes),                     
                               uy=NA, ly=NA, ul=NA, ll=NA, data=NA)
                         
            }
        } else {
            mes <- "(no estimated, no convergence model)"
            ans[[i]] <- list(lloq=NA, uloq=NA, 
                    method=paste0("coefficient variation <=",
                            max.cv, mes),                     
                            uy=NA, ly=NA, ul=NA, ll=NA, data=NA)
        
        }
    }
    class(ans) <- c("loq_cv", "loq")  
    return(ans)
}

