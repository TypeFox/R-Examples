##########################
##########################
##make pems
##########################
##########################

#kr

#description
##########################
#functions for making pems object


#includes 
##########################
#pems (nee makePEMS)
#is.pems (nee isPEMS)
#pems.element (nee makePEMSElement)
#as.pems...


#to do
##########################

#comments
##########################





##########################
##########################
##pems   nee makePEMS
##########################
##########################

#kr 18/09/2015 ver 0.0.2

#what it does
##########################
#make pems objects from parts

#to do
##########################
#is.null handling for args
#is.wrong handling for args
#defaults for constants
##########################
#

#comments
##########################
#widely used. 
#think carefully before changing name or argument ordering


pems <- function(x, units = NULL, constants = NULL,  
                     history = NULL, ...){

#################
#currently assuming 
# x = data.frame
#################

##################
#testing
#supply a pems/return it
#might want to unpack and repack?
#################

                 if(is(x)[1]=="pems") return(x)

##################
#testing 
#allow x = vector or pems.element
##################

    

    if(is.null(units) && "units" %in% names(attributes(x)))
        units <- attr(x, "units")
    x <- as.data.frame(x)


#reported issue if data.frame[1,n]


    #data has dimension to work with
    if(!is.null(ncol(x)) && ncol(x)>0){

        #units
        if(is.null(units))
            units <- rep(NA, ncol(x))
        if(!is.data.frame(units)){
            units <- as.data.frame(t(units), stringsAsFactors = FALSE)
            units <- if(ncol(units)<ncol(x))
                         cbind(units, as.data.frame(t(rep(NA, ncol(x)-ncol(units))), stringsAsFactors = FALSE)) else
                             units[1:ncol(x)] 
            names(units) <- c(names(x), names(units), rep(NA, ncol(x)))[1:ncol(x)]
        }

    }

#to do
####################
#update constants

    if(is.null(history)) history <- list()
    extra.args <- list(...)

    #update silently?
    test <- if("silent" %in% names(extra.args))
                extra.args$silent else FALSE
    extra.args <- extra.args[names(extra.args)!="silent"]

    #history
    history <- if(test)
                   history else c(history, match.call())    

    #output
    output <- list(data = x, units = units, constants = constants, 
                   history = history)

    #add in ... args
    temp <- extra.args
    output[names(temp)] <- temp

    class(output) <- "pems"

#do I want this to be invisible?

    invisible(output)
}


makePEMS <- function(...) pems(...)








##########################
##########################
##is.pems nee isPEMS
##########################
##########################

#kr 18/09/2015 v 0.0.2

#what it does
##########################
#is.pems -two level tester

#to do
##########################
#make test more robust?

#comments
##########################
#widely used. 
#think carefully before changing name or argument ordering


is.pems <- function(x, full.test = TRUE, ...){

   #standard test
   output <- if(is(x)[1]=="pems") TRUE else FALSE
   #full.test
   if(full.test){
       if(is.null(x)) comment(output) <- "NULL" else 
           if(is(x)[1]=="pems") comment(output) <- "pems" else
               if(is.data.frame(x)) comment(output) <- "data.frame" else
                    comment(output) <- "other"
   }
   #output
   output
}

isPEMS <- function(...) is.pems(...)



########################
########################
##pems.element    nee makePEMSElement
########################
########################

pems.element <- function(x, name=NULL, units=NULL, ...){

    attr(x, "class") <- unique(c("pems.element", attr(x, "class")))
    attr(x, "name") <- name
    attr(x, "units") <- units

    invisible(x)

}

makePEMSElement <- function(...) pems.element(...)




#######################
#######################
##as.pems....
#######################
#######################

##as.pems @S3 setup 
as.pems <- function(x,...)
                  UseMethod("as.pems")

##as.pems @S3 default
as.pems.default <- function(x,...){

#might need to think about this
    if(class(x)[1]=="pems") return(x)

    stop("no 'as.pems...' method for class ", 
        class(x), call. = FALSE)

}

as.pems.data.frame <- function(x,...) pems(x,...)



