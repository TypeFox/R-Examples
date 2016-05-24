#' Estimate ICC of a data.frame
#' 
#' @description Given a \code{data.frame} in long format 
#' estimate ICC from the \code{irr} package 
#' 
#' @usage
#' intra_icc(x, id.var = c("sample", "analyte"), value.var = "dil.fitted.conc",
#'     by = NULL, ...)
#' 
#' @param x a \code{data.frame} in long format with analyte, 
#' control number and concentration variables
#' @param id.var one or more variables that identify each one of 
#' replicates samples
#' @param value.var character vector with the name of the variable 
#' to estimate icc
#' @param by character vector of the variable to stratify the icc results
#' @param ... arguments for the \code{icc} function from the \code{irr} package
#' 
#' @details 
#' The \code{\link[irr]{icc}} function is the one from the \code{irr} package.
#'
#' @seealso irr
#'
#' @return
#' A list with three objects is returned:
#' \itemize{
#' \item{\code{icc.df}}{, trasnformed \code{data.frame} in 
#' order to perform the icc analysis} 
#' \item{\code{icc.mod}}{, the model output of the icc}  
#' \item{\code{icc.value}}{, the icc value of the icc model}
#' }
#' @import irr
#' 
#' @examples
#' # Generate data.frame
#' set.seed(123)
#' controls <- sort(paste("Control", rep(1:3,4),sep=""))
#' values <-  sort(unlist(lapply(1:12, function(x)runif(1,10+x,13+x))))
#' plateno <- rep(c("plate1","plate2"),6)
#' df <- data.frame(controls,values, plateno)
#' df <- df[order(df$plateno),]
#' 
#' # Estimate ICC
#' intra_icc(df, id.var = c("controls","plateno"), 
#' value.var = "values", type="agreement",model="twoway", unit="single")
#' intra_icc(df, id.var = c("controls","plateno"), 
#' value.var = "values", by = "plateno", type="agreement",model="twoway", 
#' unit="single") 
#' 
#' 
#' @export
intra_icc <- function(x, id.var = c("sample", "analyte"), 
        value.var = "dil.fitted.conc", by = NULL, ...){

    if(!inherits(x, "data.frame")) stop ("'x' must be a data.frame")
    if(".id"%in%x[,c(id.var, value.var)]) stop(" .id variable must be changed")
    if(!is.null(by) & length(by)>1) stop(" 'by' must specify only one variable")

    x$.id <- apply(x[, id.var], 1, function(y) paste(y, collapse="_"))  
    selec_vars <- unique(c(".id" , id.var, by))
    x <- x[, c(selec_vars, value.var)]
    def_info <- unique(x[,selec_vars])
    def_cont <- unique(x[,".id"])

    newx <- lapply(def_cont, function(y) c(y, x[x$.id==y, value.var]))
    newx <- ldply(newx, "rbind.fill")

    newvars <- paste(value.var, 1:(ncol(newx) -1),sep="_" )
    names(newx) <- c(".id",  newvars)
    newx <- merge(newx, def_info, by = ".id", all.x=TRUE, sort=FALSE)
    newx <- newx[,-1]
    for(i in newvars) newx[,i] <- as.numeric(newx[,i])

    if(!is.null(by)){
        by_cat <- unique(newx[,by])
        icc_model <- lapply(by_cat, 
                            function(x) icc(newx[newx[,by]==x,newvars], ...))
        names(icc_model) <- by_cat    
        icc_value <- lapply(icc_model, function(x) x$value)
        names(icc_value) <- by_cat    
    } else {
        icc_model <- icc(newx[, newvars], ...)
        icc_value <- icc_model$value
    }
    ans <- list(icc.df = newx, icc.mod = icc_model, icc.value = icc_value)
    return(ans)
}


