#' Extract a column from an event history object.
#' 
#' Extract a column from an event history object, as obtained with the function
#' \code{\link{Hist}}.
#' 
#' Since objects of class \code{"Hist"} are also matrices, all columns are
#' numeric or integer valued. To extract a correctly labeled version, the
#' attribute \code{states} of the object is used to generate factor levels.
#' 
#' @aliases getEvent 
#' @param object Object of class \code{"Hist"}.
#' @param mode Return mode. One of \code{"numeric"}, \code{"character"}, or
#' \code{"factor"}.
#' @param column Name of the column to extract from the object.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{Hist}}
#' @keywords survival
#' @examples
#' 
#'   dat= data.frame(time=1:5,event=letters[1:5])
#'   x=with(dat,Hist(time,event))
#'   ## inside integer
#'   unclass(x)
#'   ## extract event (the extra level "unknown" is for censored data)
#'   getEvent(x)
#'
#' @export
getEvent <- function(object,mode="factor",column="event"){
    model <- attr(object,"model")
    if (model=="multi.state")
        stop("Dont know how to extract events from a multi.state model")
    ## cens.code <- attr(object,"cens.code")
    states <- attr(object,"states")
    if (match(column,colnames(object),nomatch=0)==0){
        warning("Object '", class(object),"' does not have this element: ",column,". Returning NULL.")
        return(NULL)
    }
    else{
        E <- factor(as.vector(object[,column]),
                    levels=1:(length(states)+1),
                    labels=c(as.character(states),"unknown"))
        switch(mode,"character"=as.character(E),"numeric"=as.numeric(E),E)
    }
}
