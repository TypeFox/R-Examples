########################
########################
##check... functions
########################
########################

#in place
#################
#checkInput
#checkOption
#checkPEMS
#checkOuput
#checkIfMissing


#TO DO
################
#if kept 
#extend checkIfMissing to other check... arguments
#except checkOption?
################
#rewrite/tidy
#checkInput, checkOption, checkPEMS
################
#check null setting on ...Input fun.name
#could this be "checkInput" and drop some code? 
################
#

#check out
#replacing ans with rarer term outside checkInput

#questions
###################
#should there be a pems/data.name in checkPEMS?
#make the get name bit of checkInput a dedicated functions to do this?
#





########################
########################
##checkInput
########################
########################

#version 0.2.0
#karl 17/09/2010

#this can be simplified once 
#pems and pems.element full set up


checkInput <- function(input = NULL, data = NULL, 
                       input.name = NULL, fun.name = NULL, 
                       if.missing = c("stop", "warning", "return"),
                       output = c("input", "test.result"),
                       ... ){

#######################
#check defaults present/sensible
#######################


    input.name <- if(!is.character(input.name)) "input" else input.name[1]

    fun.name <- if(!is.character(fun.name)) "checkInput" else fun.name[1]

    pems <- checkPEMS(data)

    output <- checkOption(output[1], eval(formals(checkInput)$output), 
                          "output", "allowed outputs", 
                          fun.name = "checkInput")

    if.missing <- checkOption(if.missing[1], eval(formals(checkInput)$if.missing), 
                              "if.missing", "allowed if.missings", 
                              fun.name = "checkInput")

#######################
#does not seem to like
#try in if loop anymore
#hence test <- ... below
#test used in both
#######################
#later stuff might now 
#redundant....
#######################
#need to look into both
#######################

    test <- try(input, silent=TRUE)
    if(is.null(test)){
         if(if.missing == "stop")
              stop(paste("\t In ", fun.name,"(...) no ", input.name, " set", sep=""), 
                  call. = FALSE, domain = NA)
         if(if.missing== "warning")
             warning(paste("In ", fun.name,"(...) no ", input.name, " set", sep=""),
                     "\n\t [returning NULL]", 
                     call. = FALSE, domain = NA)
         return(NULL)
    }

###############
#real bit 
#to supercede hijack
###############

    if(is(test)[1] != "try-error"){
        if(is.character(input) && length(input)==1 && !is.null(data)){
            test2 <- try(data[input], silent=TRUE)
            if(is(test2)[1] != "try-error") 
                if(output=="test.result") return(TRUE) else return(test2)
        }
        if(output=="test.result") return(TRUE) else return(test)
    }

    ###########################
    #sneaky chase back of input
    ###########################
    ..x.ans <- input.name
    for(i in length(sys.frames()):1){
        temp <- eval(parse(text=
                     paste("deparse(substitute(", ..x.ans, ", sys.frames()[[i]]))", sep="")
                     ))
        if(length(temp)<2 & temp[1] !=..x.ans & temp[1] !="") ..x.ans <- temp
    }

    #tidy character inputs
    ..x.ans <- gsub("\"", "", ..x.ans)

    #if pems units available add in
    #unit NA in pems?
    #rethink this anyway

    #also want name handled in same way as below?

    units <- if(isPEMS(pems))
                 as.character(units(pems)[1, ..x.ans]) else NULL
    if(length(units)<1) units <- NULL

    
    ##########################
    #check data is there
    ##########################
    #not strictly need if used conventionally
    #but makes error tracking easier 
    if(is(try(eval(pems), silent=TRUE))[1]=="try-error")
        if(output=="test.result") return(FALSE) else
            stop(paste("\t In ", fun.name,"(...) pems/data source not found", sep=""), 
                call. = FALSE, domain = NA)

    #########################
    #recover and test input
    #########################
    #first for local case
    #

    if(length(temp)>1 || is.null(data)) {
         text <- temp 
         temp <- try(eval(parse(text = text)), silent=TRUE)
    } else { 
         text <- paste("pems$", ..x.ans, sep="") 
         text <- try(eval(parse(text = text)), silent=TRUE)
         temp <- if(class(text)[1]=="try-error" || is.null(text))
                     try(eval(parse(text = temp)), silent=TRUE) else 
                     text
    }

##########################
#this needs thinking about because it puts the unit of 
#that name in the pems source on that case
#what we really need is if input comes from 
#source and units there add them
#make object pems2
##########################


#print(temp)

#print(is(data))

#return(try(eval(parse(text = temp)), silent=TRUE))

#    temp <- if(length(temp)>1)
#                try(eval(parse(text = temp)), silent=TRUE) else
#                try(eval(parse(text = paste("pems$, ", ..x.ans, 
#                                            sep = "")))
#                    , silent = TRUE)

    if(is.function(temp)){
        temp <- "temp"
        class(temp) <- "try-error" #catch functions as well!
    }



    if(output=="test.result"){
        temp <- if(is(temp)[1]=="try-error") FALSE else TRUE
        attr(temp, "name") <- ..x.ans
        if(!is.null(units))
            attr(temp, "units") <- units
        return(temp)
    }

    if(output=="input"){
        if(is(temp)[1]=="try-error")
            switch(if.missing,
                return = { return(NULL) },
                warning = { warning(paste("In ", fun.name,"(...) ", input.name, " '", ..x.ans[1], "' not found", sep=""),
                    ", but ignoring", "\n\t [returning NULL]", call. = FALSE, domain = NA)
                return(NULL) },
                stop = { stop(paste("\t In ", fun.name,"(...) ", input.name, " '", ..x.ans[1], "' not found", sep=""),
                    call. = FALSE, domain = NA) } ) else {
            if(length(temp)==1) temp <- as.vector(temp)
            attr(temp, "name") <- ..x.ans
            if(!is.null(units))
                attr(temp, "units") <- units
###################
#new bit testing
###################
#kills pems.element posixct
#
#            class(temp) <- "pems.element"
            class(temp) <- unique(c("pems.element", class(temp)))
###################
            return(temp) } 
    }

    stop(paste("\t In ", fun.name,"(...) unexpected error", sep=""),
        "\n\t [please contact pems admin]", call. = FALSE, domain = NA)  

}







########################
########################
##checkInput2
########################
########################

#side steps hijack

#version 0.3.0
#karl 17/09/2010

#this can be simplified once 
#pems and pems.element full set up


checkInput2 <- function(input = NULL, data = NULL, 
                       input.name = NULL, fun.name = NULL, 
                       if.missing = c("stop", "warning", "return"),
                       output = c("input", "test.result"),
                       ... ){


#new
    temp <- try(comment(input), silent=TRUE)
    if(class(temp)[1]!="try-error" && !is.null(temp) && temp[1]=="hijack")
         return(input)
#note
#this current ignores missing and output
#


#######################
#check defaults present/sensible
#######################


    input.name <- if(!is.character(input.name)) "input" else input.name[1]

    fun.name <- if(!is.character(fun.name)) "checkInput" else fun.name[1]

    pems <- checkPEMS(data)

    output <- checkOption(output[1], eval(formals(checkInput)$output), 
                          "output", "allowed outputs", 
                          fun.name = "checkInput")

    if.missing <- checkOption(if.missing[1], eval(formals(checkInput)$if.missing), 
                              "if.missing", "allowed if.missings", 
                              fun.name = "checkInput")
#new
suppressWarnings(
#this kills warning below
#

    if(is.null(try(input, silent=TRUE))){
         if(if.missing == "stop")
              stop(paste("\t In ", fun.name,"(...) no ", input.name, " set", sep=""), 
                  call. = FALSE, domain = NA)
         if(if.missing== "warning")
             warning(paste("In ", fun.name,"(...) no ", input.name, " set", sep=""),
                     "\n\t [returning NULL]", 
                     call. = FALSE, domain = NA)
         return(NULL)
    }

#new
)
#


    ###########################
    #sneaky chase back of input
    ###########################
    ..x.ans <- input.name

#new
suppressWarnings(
#
    for(i in length(sys.frames()):1){
        temp <- eval(parse(text=
                     paste("deparse(substitute(", ..x.ans, ", sys.frames()[[i]]))", sep="")
                     ))
        if(length(temp)<2 & temp[1] !=..x.ans & temp[1] !="") ..x.ans <- temp
    }


#new
)
#

    #tidy character inputs
    ..x.ans <- gsub("\"", "", ..x.ans)

    #if pems units available add in
    #unit NA in pems?
    #rethink this anyway

    #also want name handled in same way as below?

    units <- if(isPEMS(pems))
                 as.character(units(pems)[1, ..x.ans]) else NULL
    if(length(units)<1) units <- NULL

    
    ##########################
    #check data is there
    ##########################
    #not strictly need if used conventionally
    #but makes error tracking easier 
    if(is(try(eval(pems), silent=TRUE))[1]=="try-error")
        if(output=="test.result") return(FALSE) else
            stop(paste("\t In ", fun.name,"(...) pems/data source not found", sep=""), 
                call. = FALSE, domain = NA)

    #########################
    #recover and test input
    #########################
    #first for local case
    #

    if(length(temp)>1 || is.null(data)) {
         text <- temp 
         temp <- try(eval(parse(text = text)), silent=TRUE)
    } else { 
         text <- paste("pems$", ..x.ans, sep="") 
         text <- try(eval(parse(text = text)), silent=TRUE)
         temp <- if(class(text)[1]=="try-error" || is.null(text))
                     try(eval(parse(text = temp)), silent=TRUE) else 
                     text
    }

##########################
#this needs thinking about because it puts the unit of 
#that name in the pems source on that case
#what we really need is if input comes from 
#source and units there add them
#make object pems2
##########################


    if(is.function(temp)){
        temp <- "temp"
        class(temp) <- "try-error" #catch functions as well!
    }



    if(output=="test.result"){
        temp <- if(is(temp)[1]=="try-error") FALSE else TRUE
        attr(temp, "name") <- ..x.ans
        if(!is.null(units))
            attr(temp, "units") <- units
        return(temp)
    }

    if(output=="input"){
        if(is(temp)[1]=="try-error")
            switch(if.missing,
                return = { return(NULL) },
                warning = { warning(paste("In ", fun.name,"(...) ", input.name, " '", ..x.ans[1], "' not found", sep=""),
                    ", but ignoring", "\n\t [returning NULL]", call. = FALSE, domain = NA)
                return(NULL) },
                stop = { stop(paste("\t In ", fun.name,"(...) ", input.name, " '", ..x.ans[1], "' not found", sep=""),
                    call. = FALSE, domain = NA) } ) else {
            if(length(temp)==1) temp <- as.vector(temp)
            attr(temp, "name") <- ..x.ans
            if(!is.null(units))
                attr(temp, "units") <- units

#new
           comment(temp) <- "hijack"
#
            class(temp) <- "pems.element"
            return(temp) } 
    }

    stop(paste("\t In ", fun.name,"(...) unexpected error", sep=""),
        "\n\t [please contact pems admin]", call. = FALSE, domain = NA)  

}










########################
########################
##checkOption
########################
########################

#version 0.2.0
#karl 17/09/2010



checkOption <- function (option=NULL, allowed.options=NULL,
  option.name = "option", allowed.options.name = "allowed.options",
  partial.match=TRUE, fun.name = "checkOption", if.missing = c("stop", "warning", "return"),
  output = c("option", "test.result"), ...){

#check options
#################
#version 0.0.1
#karl 16/09/2010

#notes
#################

#to do
#################

#################
#requested test?
#################
fun.name <- if(!is.character(fun.name)) 
    eval(formals(checkOption)$fun.name)[1] else fun.name[1]
option.name <- if(!is.character(option.name)) 
    eval(formals(checkOption)$option.name)[1] else option.name[1]
allowed.options.name <- if(!is.character(allowed.options.name)) 
    eval(formals(checkOption)$allowed.options.name)[1] else allowed.options.name[1]

#################
#if missing
#################
if (is.null(option))
   stop(paste("\t In ", fun.name,"(...) no ", option.name," set", sep=""), call. = FALSE, domain = NA)
if (is.null(allowed.options))
   stop(paste("\t In ", fun.name,"(...) no ", allowed.options.name, " set", sep=""), call. = FALSE, domain = NA)
if (is.call(allowed.options)) 
   allowed.options <- eval(allowed.options)
if (is.null(if.missing))
   stop(paste("\t In checkOption(...) no ", if.missing," set", sep=""), call. = FALSE, domain = NA)


#################################
#only first cases considered for:
#################################
option <- option[1]

##################################
#local test for output
##################################
#recursive hell for output if we don't do this
foo <- function(a,b,c, d = "checkOption"){
   ans <- if(partial.match) pmatch(a,b) else match(a,b)
   if(is.na(ans))
       stop(paste("\t In ", d, "(...) set ", c, " '", a, "' not known", sep=""),
           "\n\t [suggest one of: ", paste(b, sep=", ", collapse=", "), "]", 
           call. = FALSE, domain = NA)
   eval(b)[ans]
}

#test if.missing
if.missing <- foo(if.missing[1], eval(formals(checkOption)$if.missing), "if.missing")

#test outputs
output <- foo(output[1], eval(formals(checkOption)$output), "output")


#################################
#test options
#################################
#rationalise this?
ans <- if(partial.match) pmatch(option, allowed.options) else 
           match(option, allowed.options)


#################################
#outputs
#################################
if(output=="test.result") 
   if(is.na(ans)) return(FALSE) else {
      temp <- TRUE
      comment(temp) <- ans
      return(temp)
   }
if(output=="option")
   if(is.na(ans))
      switch(if.missing, 
          return = { return(NULL) },
          warning = { warning(paste("In ", fun.name,"(...) set ", option.name, " '", option, "' not known", 
              ", but ignoring", sep=""), "\n\t [returning NULL]", call. = FALSE, domain = NA) 
              return(NULL) },
          stop = { stop(paste("\t In ", fun.name,"(...) set ", option.name, " '", option, "' not known", sep=""),
              "\n\t [suggest one of: ", paste(allowed.options, sep=", ", collapse=", "),"]", 
              call. = FALSE, domain = NA) } ) else 
       return(allowed.options[ans]) else 
    stop(paste("\t In ", fun.name,"(...) unexpected error", sep=""),
       "\n\t [please report to pems admin]", 
       call. = FALSE, domain = NA)
}







##############################
##############################
##checkPEMS
##############################
##############################

#version 0.2.1
#karl 23/10/2010



checkPEMS <- function (data = NULL, fun.name = "checkPEMS",
                  if.missing = c("return", "warning", "stop"),
                  output = c("pems", "data.frame", "test.result"),
                  ...){

#notes
#################

#to do
#################
#see below
#

###############
#set up
###############
allowed.pems <- c("pems", "data.frame")

################
#requested test?
################

fun.name <- if(!is.character(fun.name)) 
    eval(formals(checkPEMS)$fun.name)[1] else fun.name[1]

##################
#checkPEMS output
##################
#this can not be passed back
#not a parent(..)$output issue

output <- checkOption(output[1], eval(formals(checkPEMS)$output), 
                      "output", "allowed outputs", 
                      fun.name="checkPEMS")

#################
#if.missing 
################

if.missing <- checkOption(if.missing[1], eval(formals(checkPEMS)$if.missing), 
                       "if.missing", "allowed if.missing", 
                       fun.name=fun.name)

#################
#test pems
#################

data.is <- isPEMS(data, full.test = TRUE)

################
#make if data.frame
################

if(comment(data.is)=="data.frame") data <- makePEMS(data)


#################
#output
#################

#if test
if(output=="test.result") return(data.is)

#if null, other, pems or data
if(comment(data.is)=="NULL") 
   switch(if.missing,
       return = { return(NULL) },
       warning = { warning(paste("In ", fun.name,"(...) no pems set but allowing", sep=""),
           "\n\t [returning NULL]", call. = FALSE, domain = NA)
           return(NULL) },
       stop = { stop(paste("\t In ", fun.name,"(...) no pems set", sep=""),
           "\n\t [please supply object of class: ", paste(allowed.pems, sep=" or ", collapse=" or "),"]", 
           call. = FALSE, domain = NA) } ) else 
       if(comment(data.is)=="other")
           stop(paste("\t In ", fun.name,"(...) set pems not allowed class", sep=""),
               "\n\t [suggest object of class: ", paste(allowed.pems, sep=" or ", collapse=" or "),"]", 
               call. = FALSE, domain = NA) else
           if(output=="data.frame") return(invisible(data$data)) else
               if(output=="pems") return(invisible(data)) else
                   stop(paste("\t In ", fun.name,"(...) unexpected error", sep=""),
                       "\n\t [please report to pems admin]", 
                       call. = FALSE, domain = NA)
}








###########################
###########################
##checkUnits
###########################
###########################

#version 0.2.1
#karl 23/10/2010

#issues
#########################
#possible issue rethink???
#if units not known output = "input" makes error/NULL
#


checkUnits <- function (input = NULL, units = NULL, data = NULL, 
                        input.name = NULL, fun.name = "checkUnits", 
                        if.missing = c("stop", "warning", "return"),
                        output = c("special", "units", "input", "test.result"),
                        ..., unit.conversions = NULL){

    #setups

    if.missing <- checkOption(if.missing[1], formals(checkUnits)$if.missing, 
                       "if.missing", "allowed if.missing", 
                       fun.name=fun.name)

    output <- checkOption(output[1], formals(checkUnits)$output, 
                       "output", "allowed outputs", 
                       fun.name=fun.name)
    if(output=="special")
        output <- if(is.null(units))
                      "units" else "input"

    #get input

#    ans <- checkInput(input = input, data = data, input.name = input.name, 
#                      fun.name = fun.name, if.missing = if.missing, 
#                      output = "input")


ans <- input

    input.name <- if(is.character(input.name))
                       input.name[1] else attributes(ans)$name[1]

    #make outputs

    if(output=="test.result")
        if(is.null(attributes(ans)$units[1]))
            return(FALSE) else return(TRUE)

    input.units <- attributes(ans)$units[1]

    temp2 <- if(is.null(input.name))
                "input" else paste("input '", input.name, "'", sep = "")

    if(output=="units")
        if(length(input.units) < 1){
            if(if.missing=="stop")
                stop(paste("\t In ", fun.name,"(...) no units assigned to ",
                           temp2, sep=""),
                     "\n\t [suggest assigning units in call]", 
                     call. = FALSE, domain = NA)
            if(if.missing=="warning")
                warning(paste("In ", fun.name,"(...) no units assigned to ",  
                              temp2, sep=""),
                        "\n\t [assuming right!]", 
                        "\n\t [suggest assigning units and re-running if incorrect]", 
                        call. = FALSE, domain = NA)
            return(NULL) 
        } else return(as.character(input.units)) 

    #out as input
    if(output=="input"){ 

        #no from units
        if(is.null(units)){
            if(if.missing=="stop")
                stop(paste("\t In ", fun.name,"(...) no units assigned to ",
                           temp2, sep=""),
                     "\n\t [suggest assigning units in call]", 
                     call. = FALSE, domain = NA)
            if(if.missing=="warning")
                warning(paste("In ", fun.name,"(...) no units assigned to ",  
                              temp2, sep=""),
                        "\n\t [assuming right!]", 
                        "\n\t [suggest assigning units and re-running if incorrect]", 
                        call. = FALSE, domain = NA)
            return(ans)
        }

        #no to units
        if(is.null(input.units)){
            if(if.missing=="stop")
                stop(paste("\t In ", fun.name,"(...) no units assigned to convert ",
                           temp2, " from", sep=""),
                     "\n\t [suggest assigning units to input]", 
                     call. = FALSE, domain = NA)
            if(if.missing=="warning")
                warning(paste("In ", fun.name,"(...) no units assigned to convert ",  
                              temp2, " from", sep=""),
                        "\n\t [assuming right!]", 
                        "\n\t [suggest assigning units and re-running if incorrect]", 
                        call. = FALSE, domain = NA)
            return(ans)
        }

        #if units match return
        
#####################
#testing 
#####################

       #moved forward
       if(is.null(unit.conversions)){
            unit.conversions <- ref.unit.conversions
       }

        #get all aliases of input.units
        temp <- unlist(sapply(unit.conversions, function(x) if(input.units %in% x$to) x$to))
        temp <- unique(c(input.units, temp))
        if(units %in% temp)
            return(ans)

        #replacing#
        #if(units==input.units)
        #    return(ans)

        #if new units set
        #check for unique conversion

        temp <- sapply(unit.conversions, function(x) 
                                            if(units %in% x$to & input.units %in% x$from) 
                                                TRUE else FALSE)    
###################
#replace switch
#bad idea
###################

        if(length(temp[temp])<1){
             switch(if.missing,
                   return = {return(ans)},
                   warning = {warning(paste("In ", fun.name,"(...) no unit conversion for ", input.units, 
                                            " to ", units, sep=""),
                                      paste("\n\t [assuming equivalent to default (", units, ")]", sep =""),
                                      "\n\t [suggest rescaling units and re-running if incorrect]", 
                                      call. = FALSE, domain = NA)
                              return(ans)},
                   stop = {stop(paste("\t In ", fun.name,"(...) no unit conversion for ", input.units, 
                                            " to ", units, sep=""),
                                "\n\t [suggest rescaling units and re-running if incorrect]", 
                                call. = FALSE, domain = NA)}) 
        }

        if(length(temp[temp])>1){
            warning(paste("In ", fun.name,"(...) multipe conversion methods for ", input.units, 
                          " to ", units, " conversion!", sep=""),
                    paste("\n\t [applying first]", sep =""),
                    "\n\t [suggest rescaling units and re-running if incorrect]", 
                    call. = FALSE, domain = NA)    
        }              

        temp <- unit.conversions[temp][[1]]

        ans <- temp$conversion(ans)
        attr(ans, "units") <- units
        return(ans)
    } 

    #never get to here
    stop(paste("\t In ", fun.name,"(...) unexpected error", sep=""),
         "\n\t [please report to pems admin]", 
         call. = FALSE, domain = NA)

}












##############################
##############################
##checkOutput
##############################
##############################

#version 0.2.1
#karl 23/01/2012
#     10/02/2012


#this could be tidied based 
#once pems and pems.element fully structured


checkOutput <- function (input = NULL, data = NULL,
                  input.name = NULL, fun.name = "checkOutput", 
                  if.missing = c("stop", "warning", "return"),
                  output = c("pems", "data.frame", "input", "test.result"),
                  overwrite = FALSE, ...){

#notes
#################

#to do
#################
#see below
#


    #setups

    #if.missing handling for if.missing and output?

    if.missing <- checkOption(if.missing[1], formals(checkOutput)$if.missing, 
                       "if.missing", "allowed if.missing", 
                       fun.name=fun.name)

    output <- checkOption(output[1], formals(checkOutput)$output, 
                       "output", "allowed outputs", 
                       fun.name=fun.name)

    #get input
#    ans <- checkInput(input = input, data = data, input.name = input.name, 
#                      fun.name = fun.name, if.missing = if.missing, 
#                      output = "input")

ans <- input


    input.name <- if(is.character(input.name))
                       input.name[1] else attributes(ans)$name[1]

    #standardise pems
    #only NULL or pems from here
    pems <- checkPEMS(data, fun.name = fun.name, if.missing = "return", 
                      output = "pems")

    temp <- if(is.null(input.name))
                "input" else paste ("input '", input.name, "'", sep = "")

    #outputs

    if(output == "test.results"){
        if(is.null(ans))
            return(FALSE) else return(TRUE)
    }

    #missing stop and input null stop here regardless
    if(if.missing == "stop" & is.null(ans))
        stop(paste("\t In ", fun.name,"(...) ", temp, " not supplied or NULL", sep=""),
             paste("\n\t [suggest assigning ", temp, "]", sep=""), 
             call. = FALSE, domain = NA)

    #if output input
    if(output == "input"){
        if(is.null(ans)){
            if(if.missing == "warning")
                warning(paste("In ", fun.name,"(...) ", temp, " not supplied or NULL", sep=""),
                        "\n\t [returning NULL]",
                        paste("\n\t [suggest assigning ", temp, " if required]", sep=""), 
                        call. = FALSE, domain = NA)
         }
         return(ans)
    }
    
    #if pems null make skeleton 
    temp.pems <- if(is.null(pems)) 
                     makePEMS(data.frame(NULL)) else pems

    #if ans null 
    if(is.null(ans)){
        if(if.missing == "warning"){
            if(is.null(pems)){
                warning(paste("In ", fun.name,"(...) ", temp, " and data source not supplied or NULL", sep=""),
                        paste("\n\t [returning empty ", output," ]", sep=""),
                        paste("\n\t [suggest assigning (at least) ", temp, " if required]", sep=""), 
                        call. = FALSE, domain = NA)
             } else {
                warning(paste("In ", fun.name,"(...) ", temp, " not supplied or NULL", sep=""),
                              "\n\t [returning unmodified ", output, "]",
                        paste("\n\t [suggest assigning ", temp, " if required]", sep=""), 
                        call. = FALSE, domain = NA)
             }
        }
        if(output == "data.frame") return(pemsData(temp.pems))
        if(output == "pems") return(temp.pems)
    }

    #here input is not null
    #make pems from input if pems null
    if(is.null(pems)){
        temp.pems <- makePEMS(data.frame(input = ans))

#replace with input.name


#temp
old.class <- class(temp.pems)
class(temp.pems) <- "not.pems"

        if(!is.null(attributes(ans)$name))
            names(temp.pems$data) <- attributes(ans)$name
        if(!is.null(attributes(ans)$units)){
            temp.pems$units <- data.frame(temp = attributes(ans)$units)
            names(temp.pems$units) = names(temp.pems$data)
        }

#temp
class(temp.pems) <- old.class

        if(output == "data.frame") return(temp.pems$data)
        if(output == "pems") return(temp.pems)
    } else {

        #match pems and ans

#temp
old.class <- class(pems)
class(pems) <- "not.pems"

        if(length(ans) != nrow(pems$data)){
            if(length(ans) < nrow(pems$data))
                ans <- c(ans, rep(NA, nrow(pems$data)-length(ans))) else
                pems$data[(nrow(pems$data)+1):length(ans),]<-NA

        } 


#temp
class(pems) <- old.class


#################
#add new name element created warning?
#silent if return?
#################

#temp
old.class <- class(pems)
class(pems) <- "not.pems"

        if(overwrite){
            pems$data[input.name] <- ans
            pems$units[input.name] <- if(!is.null(attributes(ans)$units))
                                          attributes(ans)$units else NA          
        } else {
            temp <- data.frame(input = ans)
            names(temp) <- input.name
            pems$data <- cbind(pems$data, temp)
            names(pems$data) <- make.unique(names(pems$data))

            temp <- data.frame(input = if(!is.null(attributes(ans)$units))
                                          attributes(ans)$units else NA )
            names(temp) <- names(pems$data)[ncol(pems$data)]
            pems$units <- cbind(pems$units, temp)
            names(pems$units) <- make.unique(names(pems$units))
        }

#temp
class(pems) <- old.class

        if(output == "data.frame") return(pemsData(pems))
        if(output == "pems") return(pems)
    }

    #shouldn't happen
    #never get to here
    stop(paste("\t In ", fun.name,"(...) unexpected error", sep=""),
         "\n\t [please report to pems admin]", 
         call. = FALSE, domain = NA)

}















##############################
##############################
##checkifMissing
##############################
##############################

#version 0.2.1
#karl 01/02/2012
#     12/02/2012


#this is recent experiment 
#could be extended to package


checkIfMissing <- function(..., if.missing = c("stop", "warning", "return"), reply = NULL, 
                           suggest = NULL, if.warning = NULL, fun.name = NULL){

    #fun.name 
    if(is.null(fun.name))
        fun.name <- "checkIfMissing"

    #if.missing
    if.missing <- checkOption(if.missing[1], formals(checkOutput)$if.missing, 
                       "if.missing", "allowed if.missing", 
                       fun.name=fun.name)

    #extra.args
    extra.args <- list(...)

    #reply
    if(is.null(reply)){
        temp <- names(extra.args)
        if(is.null(reply)){
            reply <- if(is.null(temp))
                         "expected input not found" else 
                              paste("one of expected inputs ('", paste(temp, collapse="', '", sep=""),
                                    "') not found", sep="")
         }
    }

    #suggest 
    if(!is.null(suggest)){
        reply <- paste(reply, "\n\t [suggest ", suggest, "]", sep="")

    }

    case <- length(list(...))<1 || any(sapply(list(...), is.null))

    if(case){
    
        if(if.missing == "stop")
            stop(paste("\t In ", fun.name,"(...) ", reply, sep=""), 
                 call. = FALSE, domain = NA)

        if(if.missing == "warning"){
            if(!is.null(if.warning))
                reply <- paste(reply, "\n\t [", if.warning, "]", sep="")
            warning(paste("In ", fun.name,"(...) ", reply, sep=""), 
                 call. = FALSE, domain = NA)
        }

    }


}

