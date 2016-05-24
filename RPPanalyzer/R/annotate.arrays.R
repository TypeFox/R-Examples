#' Annotates columns of RPPA expression and background matrix
#'
#' The columns of RPPA expression and background matrix are annotated
#' using the rows "pad", "slide", "incubation_run" and "spotting_run" from the
#' array description
#' @param gprData containing RPPA data
#' @param slideDesc data frame with array description meta data
#' @note for internal use only
#' @author Heiko Mannsperger <h.mannsperger@dkfz.de>
#' @return


`annotate.arrays` <-
function (gprData,slideDesc){

    ## gprData = list with four matixes from read.gpr() and if necessary followed by sub.ID()
    ## slideDesc = data frame returned from read.slidedescription

    if(!all(names(gprData) %in% c("expression", "background", "Flags", "localization"))) {
        stop("The gpr data don't contain all necessary information.")
    }

    array.id <- array.id(slideDesc)
    array.table <- t(slideDesc)
    array.table <- rbind(array.table,array.id)
    colnames(array.table) <- array.id
    checkarrays <- colnames(array.table)==colnames(gprData[[1]])
    if ( all(checkarrays)){
        x <- list(expression=gprData$expression,background=gprData$background,arraydescription=array.table,Flags=gprData$Flags, localization=gprData$localization)
        return (x)
    }
    else{
        print("can not annotate expression matrix, please check slidedescription file")
    }
}

