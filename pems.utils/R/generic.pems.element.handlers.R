##########################
##########################
##generic pems.element.handlers
##########################
##########################

#kr 01/02/2012 v 0.3.1

#includes 
#(functions/code below) 
##########################
#as.data.frame.pems.element
#print.pems.element
#plot.pems.element
#units.pems.element
#summary.pems.element
#<-.units.pems.element
#[.pems.element
#[<-.pems.element



#to do
#########################
#export [<-.pems.element



as.data.frame.pems.element <- function(x, ...){

#####################
#translation to data.frame changed
#    temp <- data.frame(as.vector(x))
#####################
    temp <- as.data.frame(as.vector(x), drop=FALSE)
    names(temp) <- if(is.null(attr(x, "name")))
                       "x" else 
                           as.character(attr(x, "name")[1])
    names(temp) <- make.names(names(temp))
#################
#this attribute handling needs to be checked
#################
    attributes(temp[,names(temp)]) <- attributes(x)
    temp    

}


print.pems.element <- function (x, ...){


#to do
###################
#tidy code
#

#to think about
###################
#do print and cat statements need to be merged
#as single output?
###################
#

    #data 
    #print with default method and attributes stripped
    ans <- x

#############
#test
#############
#    class(ans) <- "default"
#    attributes(ans) <- NULL
#    print.default(ans, ...)

####################
#update to class handling for factors etc
#might be able to replace with an inherits???
#if I understood them
#    class(ans) <- class(ans)[class(ans)!="pems.element"]
####################
    if(length(class(ans))>1) class(ans) <- class(ans)[-1] else
         class(ans)[1] <- if("levels" %in% names(attributes(ans)))
                                "factor" else mode(ans)
    attributes(ans) <- attributes(ans)[!names(attributes(ans))%in% c("name", "units")]

#allows element to print as prior class
#############

    print(ans, ...)


#see note below in plot.pems.element about
#units = "[]" plotting

    #attr
    #local report
    temp2 <- if(is.null(attributes(x)$name))
                 " [unnamed]" else paste(" ", attributes(x)$name, sep = "")
    #old line
    #         temp2 <- paste(temp2, " [", attributes(x)$units, "]", sep = "")
    if(!is.null(attributes(x)$units)){
         temp <- paste(" [", attributes(x)$units, "]", sep="")
         if(temp != " []") temp2 <- paste(temp2, temp, sep = "")
    }


    temp2 <- paste(temp2, " [n = ", length(x), "]", sep = "")
    
    cat("pems.element;", temp2, "\n", sep="")

}





##################
#plot.pems.element
##################

plot.pems.element <- function (x, y = NULL, xlab = NULL, ylab = NULL, ...){


#to this is lattice?

#rethink order
#think about error handling for mismatching cases
#bad plots, etc
#output styles?

    #x reset
#############
#previous
#    class(x) <- "default"
#############
    if(length(class(x))>1) class(x) <- class(x)[-1] else
        class(x)[1] <- if("levels" %in% names(attributes(x)))
                             "factor" else mode(x)

##could not use && attributes(x)$units!=""
##in condition term for adding [unit] to labs
##(to stop []) in plots
##must be a better way of doing this
##current is retrospective

    #get x name
    if(is.null(y)){ 
        if(is.null(ylab)){
            ylab <- if(is.null(attributes(x)$name))
                        deparse(substitute(x)) else 
                            attributes(x)$name
            if(!is.null(attributes(x)$units)){
                temp <- paste(" [", attributes(x)$units, "]", sep="")
                if(!temp %in% c(" []", " [NA]")) ylab <- paste(ylab, temp, sep = "")
            }
        }
    } else {
        if(is.null(xlab)){
            xlab <- if(is.null(attributes(x)$name))
                        deparse(substitute(x)) else 
                            attributes(x)$name
            if(!is.null(attributes(x)$units)){
                temp <- paste(" [", attributes(x)$units, "]", sep="")
                if(!temp %in% c(" []", " [NA]"))xlab <- paste(xlab, temp, sep = "")
            }
        }
        if(is.null(ylab)){
            ylab <- if(is.null(attributes(y)$name))
                        deparse(substitute(y)) else 
                            attributes(y)$name
            if(!is.null(attributes(y)$units)){
                temp <- paste(" [", attributes(y)$units, "]", sep="")
                if(!temp %in% c(" []", " [NA]")) ylab <- paste(ylab, temp, sep = "")
            }
        }
    }

    plot(x = x, y = y, xlab = xlab, ylab = ylab,...)    

}




##################
#units.pems.element
##################

#need to think about this some more
units.pems.element <- function(x) attr(x, "units")

`units<-.pems.element` <- function(x, value) { 

    #could add padding to provent bad inserts being tried

    attr(x, "units") <- value 
    x
}


####################
#summary.pems.element
####################


#summary 


summary.pems.element <- function(object, ...){

    attr(object, "class") <- attr(object, "class")[attr(object, "class") != "pems.element"]
    summary(object)

}









##########################
##########################
##[.pems.element
##########################
##########################

#kr 31/04/2014 v 0.2.1

#what it does
##########################
#handles pems.element[] calls 
#etc
#

#to do
##########################
#tidy
#think about force, simplify

`[.pems.element` <- function(x, i, ..., force=TRUE, wrap=FALSE){

    #pems.element handling
    #x[1] element 1 of x, etc 
    
    #output 
    #x[i] with pems.element attributes retained

#    att <- attributes(x)
#    class(x) <- class(x)[class(x)!="pems.element"]

############
#new
    att <- attributes(x)
    old.class <- class(x)
    class(x)[1] <- "not.pems.element"
############

    if(!force){
        i <- if(is.character(i)) 
                 i[i %in% names(x)] else i[i %in% 1:length(x)]
    }
    if(wrap && is.numeric(i)){
        if(length(x) < max(i, na.rm=TRUE)) x <- rep(x, length.out=max(i, na.rm=TRUE))
        
    }
     
    x <- try(x[i], silent = TRUE)
    if(is(x)[1] == "try-error") 
      stop("In pems.element[i] 'i' unknown/unfound", 
          call. = FALSE)
    attributes(x) <- att
################
#new
    class(x) <- old.class
################    
    x

}






##########################
##########################
##[<-.pems.element
##########################
##########################

#kr 31/04/2014 v 0.2.1

#what it does
##########################
#handles pems.element[] <- calls 
#etc
#

#to do
##########################
#tidy
#think about force, simplify

`[<-.pems.element` <- function(x, i, ..., force=TRUE, wrap=FALSE, value){

    #pems.element handling
    #x[1] <- 2 replace first case of x with 2, etc 
    #output 
    #x with pems.element attributes retained and requested insert
    #if allowed

#####################
#think about this bit here
#and in above function
#####################

#    att <- attributes(x)
#    class(x) <- class(x)[class(x)!="pems.element"]

############
#new
    att <- attributes(x)
    old.class <- class(x)
    class(x) <- if(length(class(x))<2)
                    mode(x) else class(x)[-1]
############

    if(!force){
        if(length(i) != length(value))
                  stop("In pems.element[i]<-value: i/value dimensions mismatch", 
                        call. = FALSE)
    } 
     

    test <- try(x[i]<- value, silent = TRUE)
    if(is(test)[1]=="try-error")
       stop("In pems.element[i]<- value: cannot coerce value into x[i]", 
          call. = FALSE)
#    x <- test

###############
#we currently generate a warning message here
#if attributes don't match but be coerced...
###############

    #attributes(x) <- att
################
#new
    class(x)[1] <- "pems.element"
################   
    x

}














#############################
#working
#############################


#need to look at plot and print etc.
#Warning message:
#In class(ans) <- class(ans)[class(ans) != "pems.element"] :
#  NAs introduced by coercion
#compare with other pems.element[] and pems.element[]<-



