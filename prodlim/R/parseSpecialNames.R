##' Extract from a vector of character strings the names of special functions and auxiliary arguments
##'
##' Signals an error if an element has more arguments than specified by argument arguments.
##' @title Parse special terms 
##' @param x Vector of character strings.
##' @param special A character string: the name of the special argument.
##' @param arguments A vector which contains the arguments of the special function
##' @return A named list of parsed arguments. The names of the list are the special variable names, the elements
##' are lists of arguments.
##' @seealso model.design
##' @examples
##' 
##' ## ignore arguments
##' parseSpecialNames("treat(Z)",special="treat")
##' ## set default to 0
##' parseSpecialNames(c("log(Z)","a","log(B)"),special="log",arguments=list("base"=0))
##' ## set default to 0
##' parseSpecialNames(c("log(Z,3)","a","log(B,base=1)"),special="log",arguments=list("base"=0))
##' ## different combinations of order and names
##' parseSpecialNames(c("log(Z,3)","a","log(B,1)"),
##'                   special="log",
##'                   arguments=list("base"=0))
##' parseSpecialNames(c("log(Z,1,3)","a","log(B,u=3)"),
##'                   special="log",
##'                   arguments=list("base"=0,"u"=1))
##' parseSpecialNames(c("log(Z,u=1,base=3)","a","log(B,u=3)"),
##'                   special="log",
##'                   arguments=list("base"=0,"u"=1))
##' parseSpecialNames(c("log(Z,u=1,base=3)","a","log(B,base=8,u=3)"),
##'                   special="log",
##'                   arguments=list("base"=0,"u"=1))
##' parseSpecialNames("treat(Z,u=2)",
##'                   special="treat",
##'                   arguments=list("u"=1,"k"=1))
##' parseSpecialNames(c("treat(Z,1,u=2)","treat(B,u=2,k=3)"),
##'                   special="treat",
##'                   arguments=list("u"=NA,"k"=NULL))
##' ## does not work to set default to NULL:
##' parseSpecialNames(c("treat(Z,1,u=2)","treat(B,u=2)"),
##'                   special="treat",
##'                   arguments=list("u"=NA,"k"=NULL))
 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @export
parseSpecialNames <- function(x,special,arguments){
    if (missing(arguments)) {
        argnames <- NULL
    } else {
        argnames <- names(arguments)
    }
    ## it would be possible to vectorize the function with the regexp:
    ## paste("(",paste(special,collapse="|"),")\\(|)$",sep="")
    ## but this causes some
    ## confusion and extra work
    specialRegexp <- paste("^",special,"\\(|)$",sep="")
    posSpecial <- grep(specialRegexp,x,value=FALSE)
    if (length(posSpecial)>0){
        specialTerms <- strsplit(x[posSpecial],specialRegexp)
        ## if length is 1 then term is unspecial
        ## isSpecial <- sapply(listTerms,length)
        # check for further arguments
        termsWithArguments <- unlist(lapply(specialTerms,function(x){
            if (length(x)<2) NULL
            else strsplit(x[[2]],"[ ]*,[ ]*")}), recursive=FALSE)
        varnames <- lapply(termsWithArguments,function(x){x[[1]]})
        ## attr(varnames,"special.position") <- posSpecial
        ## only fish arguments if this is desired
        if (is.null(argnames)){
            out <- vector(mode="list",length(varnames))
            names(out) <- varnames
            return(out)
        }else{
            varnames <- unlist(varnames)
            if (length(problem <- grep("=",varnames,value=TRUE))>0)
                stop(paste("Problematic variable name '",problem,"'. Variable names used in special may not contain '='.",sep=""))
            givenArguments <- lapply(termsWithArguments,function(x){
                if (length(x)==1) NULL else x[2:length(x)]
            })
            names(givenArguments) <- varnames
            # {{{  parse arguments
            specialArgumentList <- lapply(givenArguments,function(args){
                if (!is.null(args)){
                    fullvalue <- strsplit(args,"=")
                    fullvalue <- lapply(fullvalue,function(x){ ## remove whitespace
                        gsub(" ","",x)
                    })
                    givennames <- sapply(fullvalue,function(x){
                        if (length(x)==1)
                            ""
                        else
                            x[[1]]
                    })
                    values <- lapply(fullvalue,function(x){
                        if (length(x)==1)
                            x[[1]]
                        else
                            x[[2]]
                    })
                    if(length(argnames)<length(args)) stop("Too many arguments for special function '",special,"'.")
                    realnames <- givennames[givennames!=""]
                    thismatch <- match(realnames,argnames,nomatch=0)
                    if (length(realnames)>0)
                        if (!all(thismatch))
                            stop("Argument(s) '",
                                 paste(realnames,collapse=", "),
                                 "' is not an argument of '",
                                 special,
                                 "'. Valid argument(s): '",
                                 paste(argnames,collapse=", "),"'.")
                    names(values) <- givennames
                    nadd <- length(argnames)-length(values)
                    if (nadd>0){
                        values <- c(values,rep(NA,nadd))
                    }
                    thatmatch <- match(argnames,names(values),nomatch=0)
                    names(values)[names(values)==""] <- argnames[thatmatch==0]
                    values <- values[argnames]
                    ## set defaults
                    values[is.na(values)] <- unlist(arguments)[is.na(values)]
                    values
                } else {
                    ## use defaults
                    arguments
                }
            })
            # }}}
            names(specialArgumentList) <- names(givenArguments)
            ## attr(specialArgumentList,"special.position") <- posSpecial
            specialArgumentList
        }
    } else{NULL}
}
