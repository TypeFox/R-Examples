########################
########################
##unit.handlers
########################
########################

#in place
#################
#getUnits
#setUnits
#convertUnits
#addUnitConversion
#addUnitAlias
#listUnitConversions



#TO DO
################
#urgent
#convert units to log in history?
#################
#fix for aliases in convertUnits
#################
#removeUnitConversion
#removeUnitAlias
##this needs a stop orphaning ids option?



#questions
###################
#


########################
########################
##getUnits
########################
########################

#version 0.2.0
#karl 17/09/2010


getUnits <- function(input = NULL, data = NULL, ..., 
                     if.missing = c("stop", "warning", "return"),
                     hijack = FALSE){

    #if.missing handling
    if.missing <- checkOption(if.missing[1], formals(setUnits)$if.missing, 
                              "if.missing", "allowed if.missings", 
                              fun.name = "getUnits")

    #get input, then units
    ans <- if(!hijack)
               checkInput(input = input, data = data, if.missing = if.missing, 
                          fun.name = "getUnits") else
               input
    checkUnits(ans, if.missing = if.missing, fun.name = "getUnits", hijack = hijack)

}







########################
########################
##setUnits
########################
########################

#version 0.2.0
#karl 17/09/2010


setUnits <- function(input = NULL, units = NULL, data = NULL, ..., 
                     if.missing = c("stop", "warning", "return"), 
                     output = c("input", "data.frame", "pems", "special"),
                     force = FALSE, overwrite = FALSE, hijack = FALSE){

    fun.name <- "setUnits"

    #output handling
    output <- checkOption(output[1], formals(setUnits)$output, 
                       "output", "allowed outputs", 
                       fun.name = "setUnits")
    if(output == "special"){
        output <- if(is.null(data))
                      "input" else if(comment(isPEMS(data)) == "other")
                                       "input" else comment(isPEMS(data))
    }

    #if.missing handling
    if.missing <- checkOption(if.missing[1], formals(setUnits)$if.missing, 
                              "if.missing", "allowed if.missings", 
                              fun.name = "setUnits")

    #units
    if(is.null(units)){
        if(if.missing=="stop")
            stop(paste("\t In ", fun.name,"(...) units not set/NULL", sep=""),
                 paste("\n\t [suggest setting units]", sep=""), 
                 call. = FALSE, domain = NA)
        if(if.missing=="warning" & !force)
            warning(paste("\t In ", fun.name,"(...) units not set/NULL", sep=""),
                    paste("\n\t [ignoring setUnits]", sep=""), 
                    paste("\n\t [suggest setUnits(..., force = TRUE) if you meant to delete units]", sep=""), 
                    call. = FALSE, domain = NA)
    }

    ans <- if(!hijack)
               checkInput(input = input, data = data, if.missing = if.missing, 
                          output = "input", fun.name = "setUnits") else
               input

    if(is.null(units)) units <- ""
    if(is.null(attributes(ans)$units) || force || as.character(attributes(ans)$units) == as.character(units)){
        #allow null/reset
        attr(ans, "units") <- if(units=="")
                                  NULL else units
    } else {
        if(if.missing=="stop")
            stop(paste("\t In ", fun.name,"(...) could not reset already set units", sep=""),
                 paste("\n\t [suggest using convertUnits to convert current to required units]", sep=""), 
                 paste("\n\t [or setUnits(..., force = TRUE) if reset really wanted]", sep=""), 
                 call. = FALSE, domain = NA)
        if(if.missing=="warning")
            warning(paste("\t In ", fun.name,"(...) already set units not reset", sep=""),
                    paste("\n\t [ignoring setUnits call]", sep=""), 
                    paste("\n\t [suggest using convertUnits to convert current to required units]", sep=""), 
                    paste("\n\t [or setUnits(..., force = TRUE) if reset really wanted]", sep=""), 
                    call. = FALSE, domain = NA)
    }
    if(output=="input")
        attr(ans, "class") <- unique(c("pems.element", attr(ans, "class")))

    checkOutput(input = ans, data = data, if.missing = if.missing, 
                fun.name = "setUnits", output = output, overwrite = overwrite) 
    
}









########################
########################
##convertUnits
########################
########################

#version 0.2.0
#karl 17/09/2010


convertUnits <- function(input = NULL, to = NULL, from = NULL, data = NULL, ..., 
                         if.missing = c("stop", "warning", "return"), 
                         output = c("input", "data.frame", "pems", "special"),
                         unit.conversions = NULL, force = FALSE, overwrite = FALSE,
                         hijack = FALSE){

    fun.name <- "convertUnits"

    #output handling
    output <- checkOption(output[1], formals(convertUnits)$output, 
                       "output", "allowed outputs", 
                       fun.name = fun.name)
    if(output == "special"){
        output <- if(is.null(data))
                      "input" else if(comment(isPEMS(data)) == "other")
                                       "input" else comment(isPEMS(data))
    }

    #if.missing handling
    if.missing <- checkOption(if.missing[1], formals(convertUnits)$if.missing, 
                              "if.missing", "allowed if.missings", 
                              fun.name = fun.name)

    ans <- if(!hijack)
               checkInput(input = input, data = data, if.missing = if.missing, 
                          output = "input", fun.name = fun.name) else 
               input

    #from handling
    temp <- checkUnits(ans, if.missing = "return", 
                       output = "units", fun.name = fun.name, hijack = hijack)
    if(is.null(from)){

        from <- temp 

     }else { 

#################################
#fix for later
#from could be an alias of temp
#################################

        if(!force && as.character(from) != as.character(temp)){
            if(if.missing=="stop")
                stop(paste("\t In ", fun.name,"(...) from/input unit mismatch", sep=""),
                     paste("\n\t [suggest confirming input units/conversion]", sep=""),
                     paste("\n\t [or convertUnits(..., force = TRUE) if you really want conversion forced]", sep=""),
                     call. = FALSE, domain = NA)
            if(if.missing=="warning")
                warning(paste("\t In ", fun.name,"(...) from/input unit mismatch", sep=""),
                        paste("\n\t [ignoring requested convertUnits]", sep=""),
                        paste("\n\t [suggest confirming input units/conversion]", sep=""),
                        paste("\n\t [or convertUnits(..., overwrite = TRUE) if you really want conversion forced]", sep=""),
                        call. = FALSE, domain = NA)
            from <- NULL
            to <- NULL                 

        }

    }

    #if both to and from not set
    if(is.null(from) & is.null(to)){
        if(if.missing=="stop")
            stop(paste("\t In ", fun.name,"(...) to and from not set, unknown or NULL", sep=""),
                 paste("\n\t [suggest setting both]", sep=""),
                 call. = FALSE, domain = NA) 
        if(if.missing=="warning")
            warning(paste("\t In ", fun.name,"(...) to and from not set, unknown or NULL", sep=""),
                    paste("\n\t [ignoring convertUnits request]", sep=""), 
                    paste("\n\t [suggest setting both]", sep=""), 
                    call. = FALSE, domain = NA)

    }

    #to handling
    if(is.null(to)){
        if(if.missing=="stop")
            stop(paste("\t In ", fun.name,"(...) to not set/NULL", sep=""),
                 paste("\n\t [suggest setting to]", sep=""), 
                 call. = FALSE, domain = NA)
        if(if.missing=="warning")
            warning(paste("\t In ", fun.name,"(...) to not set/NULL", sep=""),
                    paste("\n\t [ignoring setUnits]", sep=""), 
                    paste("\n\t [suggest setting to in call]", sep=""), 
                    call. = FALSE, domain = NA)
         to <- from
    }

    if(!is.null(from)){
        attributes(ans)$units <- as.character(from)
    }


    if(!is.null(to)){
        ans <- checkUnits(ans, to, unit.conversions = unit.conversions, hijack = hijack, 
                          fun.name = fun.name)
    }

    checkOutput(input = ans, data = data, if.missing = if.missing, 
                fun.name = fun.name, output = output, overwrite = overwrite) 
    
}















########################
########################
##addUnitConversion
########################
########################

#version 0.2.0
#karl 17/09/2010



addUnitConversion <- function(to = NULL, from = NULL, conversion = NULL, 
                              tag = "undocumented",
                              unit.conversions = ref.unit.conversions, ...,
                              overwrite = FALSE){

    #if not unit.conversion
    if(!is.list(unit.conversions))
        unit.conversions <- list()

    #check to, from and conversion are all there!
    if(any(sapply(list(to, from, conversion), is.null)))
        stop(paste("\t In addUnitConversion(...) need all of: to, from and conversion", sep=""),
             paste("\n\t [suggest setting all in call]", sep=""), 
             call. = FALSE, domain = NA)

    to <- as.character(to)
    from <- as.character(from)
    tag <- as.character(tag)

    if(length(to)<1 | length(from)<1)
        stop(paste("\t In addUnitConversion(...) to and/or are not viable ids", sep=""),
             paste("\n\t [suggest renaming]", sep=""), 
             call. = FALSE, domain = NA)

    if(is.numeric(conversion))
        eval(parse(text=
            paste("conversion <- function(x) x * ", conversion, sep="")
            ))

    if(!is.function(conversion))
        stop(paste("\t In addUnitConversion(...) conversion not viable as method", sep=""),
             paste("\n\t [check help ?addUnitConversion]", sep=""), 
             call. = FALSE, domain = NA)

     temp <- NULL

     if(length(unit.conversions)>0){

         temp <- sapply(unit.conversions, function(x) 
                                             if(to %in% x$to & from %in% x$from) 
                                                 TRUE else FALSE)

         if(length(temp[temp])>1){
             warning(paste("In addUnitConversion(...) multipe matching conversion methods!", sep=""),
                     paste("\n\t [corrupt unit.conversions?]", sep =""),
                     "\n\t [ignoring all but first]",
                     "\n\t [suggest checking sources]", 
                     call. = FALSE, domain = NA)    
         }   

     }

     if(is.null(temp) || !any(temp)){

         #no duplicate
         unit.conversions[[length(unit.conversions)+1]] <- list(to = to, from = from, tag = tag, conversion = conversion)

     } else {

         if(overwrite){

             unit.conversions[temp][[1]]$conversion <- conversion
             if(!is.character(unit.conversions[temp][[1]]$tag) || unit.conversions[temp][[1]]$tag == "" ||
                unit.conversions[temp][[1]]$tag == "undocuments")
                    unit.conversions[temp][[1]]$tag <- tag  

         } else {

            stop(paste("\t In addUnitConversion(...) existing conversion method encountered", sep=""),
                 paste("\n\t [suggest overwrite = TRUE if you really want to do this]", sep=""), 
                 call. = FALSE, domain = NA)

         }    

     }

    unit.conversions
}













########################
########################
##addUnitAlias
########################
########################

#version 0.2.0
#karl 17/09/2010



addUnitAlias <- function(ref = NULL, alias = NULL, 
                         unit.conversions = ref.unit.conversions, ...){

    #if not unit.conversion
    if(!is.list(unit.conversions))
         stop(paste("\t In addUnitAlias(...) no unit.conversion to reference", sep=""),
             paste("\n\t [suggest updating call/checking ?addUnitAlias]", sep=""), 
             call. = FALSE, domain = NA)

    #check ref, alias are all there!
    if(any(sapply(list(ref, alias), is.null)))
        stop(paste("\t In addUnitAlias(...) need all of: ref and alias", sep=""),
             paste("\n\t [suggest setting all in call]", sep=""), 
             call. = FALSE, domain = NA)

    ref <- as.character(ref)
    alias <- as.character(alias)

    if(length(ref)<1 | length(alias)<1)
        stop(paste("\t In addUnitAlias(...) ref and/or alias not viable ids", sep=""),
             paste("\n\t [suggest renaming]", sep=""), 
             call. = FALSE, domain = NA)

    temp <- FALSE

    for(i in 1:length(unit.conversions)){

        if(ref %in% unit.conversions[[i]]$to){
            unit.conversions[[i]]$to <- unique(c(unit.conversions[[i]]$to, alias))
            temp <- TRUE
        }

        if(ref %in% unit.conversions[[i]]$from){
            unit.conversions[[i]]$from <- unique(c(unit.conversions[[i]]$from, alias))
            temp <- TRUE
        }

    }

    if(!temp)
         warning(paste("In addUnitAlias(...) ref not found in look-up table", sep=""),
             paste("\n\t [no alias updates]", sep=""), 
             call. = FALSE, domain = NA)

    unit.conversions
}










########################
########################
##listUnitConversions
########################
########################

#version 0.2.0
#karl 17/09/2010



listUnitConversions <- function(unit.conversions = ref.unit.conversions, ...,
                                verbose = FALSE, to = NULL, from = NULL){


    #if not unit.conversion
    if(!is.list(unit.conversions))
         stop(paste("\t In listUnitConversions(...) no unit.conversion to reference", sep=""),
             paste("\n\t [suggest updating call/checking ?listUnitConversions]", sep=""), 
             call. = FALSE, domain = NA)

    #set up to, from
    to <- if(!is.null(to))
              as.character(to) else ""
    from <- if(!is.null(from))
                as.character(from) else ""

######################
#error if to, from no good?
######################

    if(to != "" | from != ""){
        temp <- sapply(unit.conversions, function(x) 
                                            if(to %in% x$to | from %in% x$from) 
                                                TRUE else FALSE)
        unit.conversions <- unit.conversions[temp]    
    }

    if(length(unit.conversions)<1)
         stop(paste("\t In listUnitConversions(...) no matched methods located", sep=""),
             paste("\n\t [no suggestion]", sep=""), 
             call. = FALSE, domain = NA)

    temp.fun <- if(verbose){
                    function(x) 
                        paste("TAG: ", paste(x$tag, sep ="", collapse =","),
                              "; FROM:", paste(x$from, sep ="", collapse =","),
                              "; TO:", paste(x$to, sep ="", collapse =","), sep="")
                } else {
                    function(x) 
                        paste(x$tag, sep ="", collapse =",")
                } 
                    
    sapply(unit.conversions, temp.fun) 

}




