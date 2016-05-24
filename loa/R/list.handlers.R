#in development code
#[TBC - NUMBER] functions 

#listHandler
#listUpdate
#listExpand
#listLoad

#NOTE: much borrowed from lattice 

##################
#possibles
##################
#possibly rethink drop.dots
#so it can drop stuff
#####
#should listExpand be listExpandVectors?
#because that is more descriptive
#####
#should listLoad have an output argument
#


##################
#to do
##################
#rethink ignore handling/descriptions
#could also have 




############################
############################
##listHandler
############################
############################

listHandler <- function(a, use = NULL, ignore = NULL, 
                        drop.dots=TRUE){
    
    if(drop.dots)
        a <- a[names(a) != "..."]
    if(!is.null(use))
        a <- a[names(a) %in% use]
    if(!is.null(ignore))
        a <- a[!names(a) %in% ignore]

    a

}



#####################
#####################
##localUpdate
#####################
#####################

#local function to update lists
#[in development]

listUpdate <- function(a, b, use = NULL, ignore = NULL,
                       use.a = use, use.b = use,
                       ignore.a = ignore, ignore.b = ignore, 
                       drop.dots = TRUE){

    a <- listHandler(a, use.a, ignore.a, drop.dots)
    b <- listHandler(b, use.b, ignore.b, drop.dots)

#################
#testing
#################

    if(is.null(a) & is.list(b)) return(b)
    if(is.list(a) & is.null(b)) return(a)

    if(length(names(b) > 0))
        a <- modifyList(a, b)
    a
}




############################
############################
##listExpand
############################
############################

listExpand <- function(a, ref = NULL, use = NULL, 
                       ignore = NULL, drop.dots = TRUE){

    a <- listHandler(a, use, ignore, drop.dots)

    if(is.null(ref))
        return(a)

    temp <- lapply(a, function(x){x <- if(is.vector(x) & !is.list(x)){
                             if(length(x) > 1 & length(x) < length(ref))
                                 rep(x, ceiling(length(ref)/length(x)))[1:length(ref)] else x
                             } else x 
                         })
    listUpdate(a, temp)
}



############################
############################
##listLoad
############################
############################

listLoad <- function(..., load = NULL){

    #######################
    #listLoad v0.1
    #######################
    #kr 21/03/2013
    #######################

    #loads list with associated prefixed terms 
    #and strips loaded cases from args

    extra.args <- list(...)
    output <- if(is.null(extra.args$output))
                   "all" else extra.args$output

#might want to restrict this to just characters
#current default output always as "all"
#current detault one list only

    if(!isGood4LOA(load)) 
        return(extra.args)

#if list is false don't load it

   if (is.logical(extra.args[[load]]) && !extra.args[[load]])
        return(extra.args) 

    load.cases <- grep(paste(load, "[.]", sep=""), names(extra.args), value=T)



    if(length(load.cases)>0){
        temp <- extra.args[load.cases]
        names(temp) <- gsub(paste(load, ".", sep=""), "", names(temp))

        extra.args <- extra.args[!names(extra.args) %in% load.cases]

####################
#rethink
####################
#load.list <- extra.args[[load]]
#

        if (is.logical(extra.args[[load]]) && extra.args[[load]]) 
            extra.args[[load]] <- list()
        if (is.list(temp))
            extra.args[[load]] <- listUpdate(temp, extra.args[[load]])
    }

    return(extra.args)

}




