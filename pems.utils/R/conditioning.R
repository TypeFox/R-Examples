##########################
##########################
##condition pems objects
##########################
##########################

#kr

#description
##########################
#functions to condition pems


#includes 
##########################
#cutBy
#


#to do
##########################

#comments
##########################
#cutBy could be intelligent
#i.e. choose a method based on ref and inputs
#see below
#also the output could be a list of 
#pems cut by the cut term




##########################
##########################
##cutBy
##########################
##########################

#kr 23/01/2012 v 0.0.6

#what it does
##########################
#wrapper for cutting data


#to do
##########################
#make test more robust?

#comments
##########################
#


cutBy <- function(ref = NULL, ..., data = NULL, cut.method = NULL, 
                  labels = NULL,
                  fun.name = "cutBy", hijack= FALSE){
  
    #setup
    this.call <- match.call()
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 
    if(!hijack)   
        ref <- checkInput(ref, data=data, fun.name = fun.name,   
                            if.missing = settings$if.missing,
                            unit.conversions = settings$unit.conversions)  

    if(is.null(cut.method)){
        #select suitable cut.method

######################
#this could be based on 
#ref and args in ...
##so e.g. if ref NULL and n or rows in ... 
##        use cutByROw 
#######################

        cut.method <- cutByRow

    }


    if(is.function(cut.method)){

        cut <- cut.method(ref = ref, data = data, output = "input",
                          ..., fun.name = "cutBy", hijack= TRUE)

        #check for attributes 
        #if not give it some

#new addition
#as part non-special behaviour...
        attr(cut, "class") <- unique(c("pems.element", attr(cut, "class")))
#
        if(is.null(attributes(cut)$name))
            attr(cut, "name") <- "cut"
        if(is.null(attributes(cut)$units))
            attr(cut, "units") <- "cut.region"

        #do any renaming
        if(!is.null(labels)){
            levels(cut) <- labels
        }
         
        return(calcPack(output = cut, data = data, settings = settings, 
                 fun.name = fun.name, this.call = this.call)) 
    }

    #not good
    checkIfMissing(if.missing = settings$if.missing, 
                   reply = "could not run cut.method!", 
                           suggest = "check ?cutBy if reason unclear", if.warning = "returning NULL", 
                           fun.name = fun.name)
    return(NULL)    
}



cutByRow <- function(ref = NULL, n = 4, rows = NULL, ..., data = NULL,  
                    fun.name = "cutByRow", hijack= FALSE){
  
    #setup
    this.call <- match.call()
    
    #run checks
    settings <- calcChecks(fun.name, ..., data = data)

    #get what there is 
    if(!hijack)   
        ref <- checkInput(ref, data=data, fun.name = fun.name,   
                            if.missing = settings$if.missing,
                            unit.conversions = settings$unit.conversions)  

    #both n and row should 
    if(!is.numeric(n) & !is.numeric(rows)){
        checkIfMissing(if.missing = settings$if.missing, 
             reply = "require at least numeric 'n' or 'rows' to cut something!", 
             suggest = "check ?cutBy", if.warning = "returning NULL", 
             fun.name = fun.name)
        return(NULL)  
    }

    

    if(is.numeric(n) & !is.numeric(rows)){
        if(n > length(ref)) n <- length(ref)
        rows <- floor(length(ref)/n[1])
        rows <- seq(1, length(ref), by = rows)
    }

    if(is.numeric(rows)){
        rows <- floor(rows)
        rows <- rows[rows>0 & rows<(length(ref)+1)]
        cut <- rep(0, length(ref))
        cut[rows] <- 1
        cut[1] <- 1
        cut <- cumsum(cut)
    }
 
    cut <- factor(cut)

    #add attributes
    attr(cut, "name") <- "cut"
    attr(cut, "units") <- "cut.region"

    #return
    return(calcPack(output = cut, data = data, settings = settings, 
                    fun.name = fun.name, this.call = this.call)) 

}








