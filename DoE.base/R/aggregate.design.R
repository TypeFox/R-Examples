aggregate.design <- function(x, ..., by=NULL, response=NULL, FUN="mean", postfix=NULL, replace=TRUE){
     if (!"design" %in% class(x)) stop("x must be of class design")
     di <- design.info(x)
     if (is.null(di$responselist) & is.null(by)) 
        stop("x must be a wide design (repeated measurements or parameter design), or by must be specified")
     if (is.null(di$responselist) & !is.null(by)) 
        return(aggregate.data.frame(x, by, FUN, ...))

    ## from here on, treatment of wide designs
     if (is.null(postfix)){
        if (is.character(FUN)) postfix <- FUN
        else postfix <- make.names(deparse(substitute(FUN)))
        }
     if (is.character(FUN) & length(FUN)>1) 
        stop("aggregate.design can only handle one function at a time")
     if (is.null(postfix)) postfix <- FUN
     if (!(is.character(postfix) & length(postfix)==1)) 
        stop("postfix must be a character string")
     FUN <- match.fun(FUN)
     if (is.null(response)) response <- names(di$responselist)
     if (!is.character(response)) stop("response must be a character vector of response names")
     if (!length(setdiff(response, colnames(di$responselist)))==0)
        stop("invalid response name(s)")
     aus <- x
     for (i in 1:length(response)){
        assign(paste(response[i],postfix,sep="."), apply(x[,di$responselist[,response[i]]],1,FUN))
        aus <- eval(parse(text=paste("add.response(aus,", paste(response[i],postfix,sep="."),", replace=replace)")))
        }
     #modified 30 Jan 2011; not useful to remove all responses
     di$response.names <- setdiff(design.info(aus)$response.names, unlist(di$responselist))
     design.info(aus) <- di
     aus
}