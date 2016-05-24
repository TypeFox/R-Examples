########################
########################
##pems.structure
########################
########################

#in place
#################
#pemsElement
#pemsData
#pemsConstants
#pemsHistory
#pemsDim



#TO DO
################

#pemsConstants
#tidy
#document
#namespace





#questions
###################
#is this better than check...
#do functions need a test that first element is pems, etc?
#




########################
########################
##pemsElement
########################
########################

#version 0.2.0
#karl 17/09/2010

pemsElement <- function(element, pems=NULL, ..., 
         fun.name = "pemsElement", if.missing = "stop",
         element.name = deparse(substitute(element))){

   #reorder this later

    if(is.null(pems)) {
       ans <- try(element, silent = TRUE)
       if(is(ans)[1] == "try-error" || is.function(ans)){ 
           checkIfMissing(if.missing = if.missing,
                          reply = paste("element '", element.name[1], "' not found", sep=""),
                          suggest = "checking call arguments", 
                          if.warning = NULL, 
                          fun.name = fun.name)
           ans <- NULL
       }
       if(!is.null(ans))
           if(is.null(attributes(ans)$name)) attr(ans, "name") <- element.name
       return(ans)   
    }

    pems <- checkPEMS(pems)
    class(pems) <- "not.pems"

    ans <- try(pems$data[,element.name], silent = TRUE)
    units <- try(pems$units[1,element.name], silent = TRUE)

    if(is(ans)[1] == "try-error" || is.function(ans)){ 
        ans <- try(element, silent = TRUE)
        if(is(ans)[1] == "try-error" || is.function(ans)){ 
             checkIfMissing(if.missing = if.missing,
                          reply = paste("element '", element.name[1], "' not found", sep=""),
                          suggest = "checking call arguments", 
                          if.warning = NULL, 
                          fun.name = fun.name)
             units <- NULL
             ans <- NULL
         }        
    }

    if(!is.null(ans))
        attr(ans, "name") <- element.name
    if(!is.null(ans) && !is.null(units))
        if(is.null(attributes(ans)$units)) attr(ans, "units") <- units
    class(ans) <- "pems.element"

    return(ans)

}






##############################
##############################
##pemsData
##############################
##############################


pemsData <- function(pems=NULL, ..., 
         fun.name = "pemsData", if.missing = "stop",
         pems.name = deparse(substitute(pems))){

    if(is.null(pems)){
          checkIfMissing(if.missing = if.missing,
               reply = paste("pems '", pems.name[1], "' not found", sep=""),
               suggest = "checking call arguments", 
               if.warning = NULL, 
               fun.name = fun.name)
    }
    class(pems) <- "not.pems"
    pems$data

}




##############################
##############################
##pemsConstants
##############################
##############################


pemsConstants <- function(pems=NULL, ..., 
         fun.name = "pemsConstants", if.missing = "stop",
         pems.name = deparse(substitute(pems))){

    if(is.null(pems)){
          checkIfMissing(if.missing = if.missing,
               reply = paste("pems '", pems.name[1], "' not found", sep=""),
               suggest = "checking call arguments", 
               if.warning = NULL, 
               fun.name = fun.name)
    }
    class(pems) <- "not.pems"
    pems$constants

}





##############################
##############################
##pemsHistory
##############################
##############################


pemsHistory <- function(pems=NULL, ..., 
         fun.name = "pemsHistory", if.missing = "stop",
         pems.name = deparse(substitute(pems))){

    if(is.null(pems)){
          checkIfMissing(if.missing = if.missing,
               reply = paste("pems '", pems.name[1], "' not found", sep=""),
               suggest = "checking call arguments", 
               if.warning = NULL, 
               fun.name = fun.name)
    }
    class(pems) <- "not.pems"
    pems$history

}








##############################
##############################
##pemsDim
##############################
##############################


pemsDim <- function(pems=NULL, ..., 
         fun.name = "pemsDim", if.missing = "stop",
         pems.name = deparse(substitute(pems))){

    if(is.null(pems)){
          checkIfMissing(if.missing = if.missing,
               reply = paste("pems '", pems.name[1], "' not found", sep=""),
               suggest = "checking call arguments", 
               if.warning = NULL, 
               fun.name = fun.name)
    }
    class(pems) <- "not.pems"
    dim(pems$data)

}
