lapplyBy <- function (formula, data = parent.frame(), FUN, keep.groupid=FALSE)
{
    ddd <- splitBy(formula, data = data)
    att <- attributes( ddd )
    grp <- unique( att$grps )
    ddd <- lapply(ddd, FUN)
    ddd <- ddd[grp] ## probably not necessary
    if (keep.groupid){
        attr(ddd, "groupid") <- att$groupid
        ##attributes( ddd ) <- att
    }
    #class( ddd ) <- c("lapplyByData", "list")
    ddd
}

## print.lapplyByData <- function(x,...){
##     oclass <- class(x)
##     class(x) <- "listof"
##     print(x)
##     class(x) <- oclass
##     return(invisible(x))
## }

## summary.lapplyByData <- function(object, ...){
##     zz <- attributes( object )
##     class(zz) <- "lapplyBySummary"
##     zz
## }

## print.lapplyBySummary <- function(x, ...){
##     print(x$groupid)
## }
