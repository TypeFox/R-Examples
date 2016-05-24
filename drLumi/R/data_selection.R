#' Extracts samples information from a lum_export object or data.frame
#'
#' @description
#' Extracts from a \code{lum_export} or \code{data.frame} object controls, 
#' dilutions points, background and samples to be calibrated. 
#' These files can be merged to an expected concentration and flag dataset.
#' 
#' @usage
#' data_selection(x, ecfile = NULL, flagsfile = NULL, backname = "Back",
#'     stanname = "Stan", posname = "Con", unkname = NULL,
#'     byvar.ecfile = c("analyte", "sample"), byvar.flagsfile = c("well",
#'     "analyte"), fsample = "sample", fanalyte = "analyte", fbatch = "plate",
#'     ...)
#' 
#' @param x a \code{lum_export} object or \code{data.frame} with 
#' samples information.
#' @param ecfile a \code{data.frame} or CSV path file with the 
#' expected concentration in order to merge to \code{x}. 
#' @param flagsfile a \code{data.frame} or CSV path with the flags 
#' information in order to merge to \code{x}. 
#' @param backname character vector or list of the background samples 
#' identification of \code{x}.
#' @param stanname character vector or list of the standard points 
#' identification of \code{x}.
#' @param posname character vector or list  of the positive controls 
#' identification of \code{x}.
#' @param unkname character vector or list of the samples identification 
#' of \code{x}. Default \code{NULL}.
#' @param byvar.ecfile character vector of the merging variable(s) 
#' of \code{x} and \code{ecfile}.
#' @param byvar.flagsfile character vector of the merging variable(s) 
#' of \code{x} and \code{flagsfile}.
#' @param fsample character vector of the name of the sample variable.
#' @param fanalyte character vector of the name of the analyte variable.
#' @param fbatch character vector of the name of the Batch variable.
#' @param ... other options. Ignored.
#' 
#' @details Default method for identifying background, standard and 
#' positives samples is to define a character vector of length one and 
#' apply \code{\link{agrep}} functions in order to extract databases. 
#' Samples to be calibrated (unknowns) are identified as the 
#' remaining ones (default \code{NULL}).
#' 
#' If the arguments are defined as a list the function will 
#' subset exactly that information from \code{x} object. 
#' 
#' The expected concentration file is merged based on \code{byvar.ecfile}. 
#' Only applies when \code{ecfile} is not \code{NULL}. Same applies for flags file. 
#' A variable named 'flag' must be in flags data in order to perform the merge.
#' 
#' @return The name of the batch in a list format with the following 
#' components: background, standard, positive and unknowns datasets.
#' 
#' @import Hmisc
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' dat <- subset(mfidata,plate=="plate_1" & analyte=="FGF")
#' 
#' # Example 1
#' sdf <- data_selection(dat)
#' 
#' lapply(sdf$plate_1, function(x) head(x))
#' 
#' # Example 2 (merge ecdata)
#' sdf <- data_selection(dat, ecfile = ecdata, 
#'              byvar.ecfile=c("analyte","sample"))
#' 
#' lapply(sdf$plate_1, function(x) head(x))
#' 
#' # Example 3 (extract specific samples names with list)
#' sdf <- data_selection(dat, 
#'              stanname=list("Standard10"), 
#'              backname = list("Background0"),
#'              posname = list("Control1","Control2"), 
#'              unkname = list("B_sid_13_CSP"))
#'                 
#' lapply(sdf$plate_1, function(x) head(x))
#' 
#' # Example 4 (extract aproximate names samples)
#' sdf <- data_selection(dat, 
#'              stanname="Standard1", 
#'              backname = "Background0", 
#'              posname = "Control1", 
#'              unkname = "B_sid_13_CSP")
#' 
#' lapply(sdf$plate_1, function(x) (x))  
#' 
#' 
#' @export
data_selection <- function(x, ecfile = NULL, flagsfile = NULL, 
                    backname = "Back", 
                    stanname = "Stan", posname = "Con", unkname = NULL, 
                    byvar.ecfile = c("analyte","sample"), 
                    byvar.flagsfile = c("well","analyte"),  
                    fsample = "sample", fanalyte = "analyte", 
                    fbatch = "plate", ...){

    if(!inherits(x, c("lum_export","data.frame"))){ 
        stop("'x' must be 'lum_export' or 'data.frame' object")
    }
    if(!is.null(ecfile) & !inherits(byvar.ecfile,"character")){
        stop("'byvar.ecfile' must be a character vector")
    }
    if(!is.null(ecfile) & !inherits(byvar.flagsfile,"character")){
        stop("'byvar.ecfile' must be a character vector")
    }
    if(!is.null(ecfile)){
        if(!inherits(ecfile, "data.frame")){
            ec <- try( read.csv(ecfile,as.is=TRUE)  , silent=TRUE)
            if(inherits(ec, "try-error")){
                warning("No possible to read 'ecfile'.")
            }  
        } else {
            ec <- ecfile
        }  
    } else{
        ec <- NA
    }
    
    if(!is.null(flagsfile)){
        if(!inherits(flagsfile, "data.frame")){
            flags <- try( read.csv(flagsfile,as.is=TRUE)  , silent=TRUE)
            if(inherits(flags, "try-error")){ 
                warning("No possible to read Flag file.")
                flags <- NA
            }  
        } else {
            flags <- flagsfile
        }  
    } else{
        flags <- NA
    }    

    if(inherits(x,"lum_export")){
        batchname <- x$name_batch
        x <- data.frame(lapply(x$well, as.character), stringsAsFactors=FALSE)
    } else{
        batchname <- as.character(unique(x[,fbatch]))
        if(length(batchname)>1) stop("Only one batch can be analyzed")
    } 
    
    if(fsample%nin%names(x)) stop("No 'fsample' in 'x' object" )
    if(fanalyte%nin%names(x)) stop("No 'fanalyte' in 'x' object" )  
    if(fbatch%nin%names(x)) stop("No 'fbatch' in 'x' object" )

    newdata <- x
    if(inherits(flags, "data.frame")){    
        newdatax <- try(merge(newdata, flags[,c(byvar.flagsfile,"flag")], 
                        by = byvar.flagsfile, 
                        all.x=TRUE, sort=FALSE), silent=TRUE)    
        if(inherits(newdatax, "try-error")){
            warning("No possible to merge 'flagsfile' with 'x' file.")
        } else {
            newdata <- newdatax
            newdata$flag <- ifelse(is.na(newdata$flag),"",newdata$flag )   
        }
    } 
    samples_names <- as.character(unique(newdata[,fsample]))

    if(!is.null(backname)){
        if(inherits(backname, "list"))  back_sel <- unlist(backname)
        if(inherits(backname, "character")){
            back_sel <- agrep(backname, samples_names, value = TRUE)  
        }  
        back_df <- newdata[newdata[,fsample]%in%back_sel,]  
        } else {
            back_sel <- NULL  
            back_df <- newdata[newdata[,fsample]%in%back_sel,]  
        }

    if(!is.null(stanname)){
        if(inherits(stanname, "list"))  stan_sel <- unlist(stanname)
        if(inherits(stanname, "character")){
            stan_sel <- agrep(stanname, samples_names, value = TRUE)    
        }
        stan_df <- newdata[newdata[,fsample]%in%stan_sel,]  
    } else {
        stan_sel <- NULL  
        stan_df <- newdata[newdata[,fsample]%in%stan_sel,]  
    }

    if(!is.null(posname)){
        if(inherits(posname, "list"))  pos_sel <- unlist(posname)
        if(inherits(posname, "character")){
            pos_sel <- agrep(posname, samples_names, value = TRUE)    
        }
        pos_df <- newdata[newdata[,fsample]%in%pos_sel,]  
    } else {
        pos_sel <- NULL  
        pos_df <- newdata[newdata[,fsample]%in%pos_sel,]  
    }  

    if(inherits(unkname, "list"))  unk_sel <- unlist(unkname)
    if(inherits(unkname, "character")){
        unk_sel <- agrep(unkname, samples_names, value = TRUE)
    }
    if(is.null(unkname)){
        snames <- samples_names
        unk_sel <- snames[snames%nin%c(pos_sel,stan_sel, back_sel)]
    }
    unk_df <- newdata[newdata[,fsample]%in%unk_sel,] 

    if(!inherits(ec, "try-error") & inherits(ec,"data.frame")){
        stan_df2 <- try(merge(stan_df, ec, 
                        by = byvar.ecfile , 
                        all.x=TRUE, sort=FALSE),silent=TRUE)
        if(inherits(stan_df2, "try-error")){
            warning("No possible to merge 'ecfile' with 'standard' file.")
        } else{
            stan_df <- stan_df2
        }  
        back_df2 <- try(merge(back_df, ec, by = byvar.ecfile, 
                    all.x=TRUE, sort=FALSE), silent=TRUE)
        if(inherits(back_df2, "try-error")){
            warning("No possible to merge 'ecfile' with 'background' file.")
        } else {
            back_df <- back_df2
        }
        pos_df2 <- try(merge(pos_df, ec, by = byvar.ecfile, 
                    all.x=TRUE, sort=FALSE), silent=TRUE)
        if(inherits(pos_df2, "try-error")){
            warning("No possible to merge 'ecfile' with 'positives' file.")
        } else{
            pos_df <- pos_df2
        }
        unk_df2 <- try(merge(unk_df, ec, by = byvar.ecfile, 
                    all.x=TRUE, sort=FALSE), silent=TRUE)
        if(inherits(unk_df2, "try-error")){
            warning("No possible to merge 'ecfile' with 'unknown' file.")
        } else{
            unk_df <- unk_df2
        }
    }
    # if(nrow(back_df)==0){
    # back_df <- data.frame(analyte = unique(newdata[,fanalyte]))
    # back_df[names(newdata)[names(newdata)%nin%fanalyte]] <- NA
    # back_df$ec <- NA
    # }
    sta_selec <- stan_df[, fanalyte]!="Total Events"
    standard <- stan_df[sta_selec,]
    pos_selec <- pos_df[, fanalyte]!="Total Events"
    positive <- pos_df[pos_selec,]
    back_selec <- back_df[, fanalyte]!="Total Events"
    background <- back_df[back_selec,]  
    unk_selec <- unk_df[, fanalyte]!="Total Events"
    unknowns <- unk_df[unk_selec,]  
    ans <- list()
    ans[[batchname]] <- list(background = background, standard = standard, 
                        positive = positive, unknowns = unknowns, 
                        name_batch = batchname)
    return(ans)
}





