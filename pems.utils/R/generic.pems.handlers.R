##########################
##########################
##generic pems handlers
##########################
##########################


#might to archive a copy of this and then tidy
#because it is a real mess at the moment


#kr 01/02/2012 v 0.3.1

#includes 
#(functions/code below) 
##########################
#as.data.frame.pems
#dim.pems
#############->nrow.pems
#############->ncol.pems
#[.pems
#[<-.pems
#[[.pems
#[[.<-pems
#$.pems
#$<-.pems
#with.pems
#subset.pems
#print.pems
#plot.pems
#names.pems
#names<-.pems
#summary.pems
#head.pems
#tail.pems
#units.pems
#units<-.pems
#

#to do
##########################
#so much
#see individual methods 
#all need work

# verbose printing? 
# but needs to be hidden when dumped 
# to object
#
# head.pems <- function(x, n =6,...) x[1:n,]
# tail

#as.list.pems <- function(x, ...) unclass(x)
#as.pems.list <- function(x, ...) do.call(makePEMS, x)
#but latter would need a method make

#


#comments
##########################
#



##########################
##########################
##as.data.frame.pems
##########################
##########################

#pull data.frame out of pems
#for lattice, lm, etc.

#think about this
#can we get a vector out rather than a data.frame?

as.data.frame.pems <- function(x, ...){

    class(x) <- "not.pems"
    x$data    

}




##########################
##########################
##dim.pems
##nrow.pems
##ncol.pems
##########################
##########################

#get pems dimensions

dim.pems <- function(x, ...) dim(as.data.frame(x))

#not needs because dim.pems gives you these
#nrow.pems <- function(x, ...) nrow(as.data.frame(x))
#ncol.pems <- function(x, ...) ncol(as.data.frame(x)) 






##########################
##########################
##[[.pems
##########################
##########################

#kr 18/08/2015 v 0.1.0

#what it does
##########################
#handles pems[[]] calls 
#for access to data, units and other tags
#

#to do
##########################
#tidy
#think about this

##need pems[[]]<- operator


`[[.pems` <- function(x, k, ...){

    #break pems
    class(x) <- "list"

    #return as list is nothing declared
    if(missing(k)) return(x)

    #select structural elements
    temp.fun <- function(k){
        ans <- try(x[[k]], silent=TRUE)
        if(class(ans)[1] == "try-error") NULL else ans
    }
    if(length(k)==1) return(temp.fun(k))
    lapply(k, temp.fun)

}




#########################
#########################
##[[<-.pems
#########################
#########################

`[[<-.pems` <- function(x, k, ..., value){

    #break pems
    old.class <- class(x)
    class(x) <- "list"

    #return as list is nothing declared
    if(missing(k) | missing(value)) return(x)

    #add in if it exists or not
    #might want to think about this?

    x[[k]] <- value 
    class(x) <- old.class
    x

}





##########################
##########################
##[.pems
##########################
##########################

#kr 06/06/2013 v 0.3.0

#what it does
##########################
#handles pems[] calls 
#etc
#

#to do
##########################
#tidy
#think about force, simplify

`[.pems` <- function(x, i, j, ..., force = FALSE, simplify = TRUE){

    ########################
    #generic pems handling
    ########################

    #x[1] element 1 as 1 col data.frame
    #x[1,] row 1
    #x[,1] element 1 
    #x[1,1] row, element 1,1 
    
    #output options
    #force=FALSE fail call if any requested elements/rows invalid
    #force=TRUE return valid element/row requests, discard invalid 
    #simplify=TRUE as pems.element if possible else pems
    #simplify=FALSE as pems 

    #force = FALSE, simplify = TRUE gives the response most like a 
    #conventional data frame
    #currently handles negs like data.frame (I think)...

    #negs are discarded if a mixed (neg and pos) vector supplied
    #and forced. Possibly a better way of handling this

    call1 <- sys.call()
    call2 <- match.call()
    check.i <- !missing(i)
    check.j <- !missing(j)

    #request pems[], pems[,], etc
    #do nothing
    if(!check.i & !check.j) return(x)

#in progress 
#current na.pad.output is functional but in testing
#might still be issues

    #force special handling
    #options omit.err.cases (default forcing)
    #        no.forcing
    #        na.pad.output
    if(is.logical(force))
        force <- if(force) "omit.err.cases" else "no.forcing"

    #make 
    old.class <- class(x)
    class(x) <- "not.pems"

    #request pems[i], etc
    #columns/elements, etc
    if(check.i && !check.j && length(as.character(call1))==length(as.character(call2))){
       j <- i
       i <- 1:nrow(x$data)
    }

    #request pems[i,], etc
    #rows, etc
    if(check.i && !check.j && length(as.character(call1))!=length(as.character(call2)))
        j <- 1:ncol(x$data)

    #request pems[,j], etc
    #columns/elements
    if(!check.i && check.j)
        i <- 1:nrow(x$data)
    
    #otherwise it is pems[i,j], etc

    #negative handling 

    #at this stage we err out on negs
    test <- character()
    if(is.numeric(j) && any(j<=0)) test <- c(test, "elements")
    if(is.numeric(i) && any(i<=0)) test <- c(test, "rows")
    if(length(test)>0)
        stop(paste("In pems[i,j]<-: neg ",
                   paste(test, sep=" and ", collapse=" and "),
                   " requested \n       [not currently allowed]",
                   sep=""), call. = FALSE)


#data.frame like neg handling

#note/known issue
#currently regards zero as negative
 
#    if(is.numeric(i) && any(i<=0))
#        if(all(i<=0))
#            i <- c(1:nrow(x$data))[!1:nrow(x$data) %in% abs(i)] else
#                 if(force=="omit.err.cases") 
#                     i <- i[i>0] else 
#                         stop("In pems[i,j]: mixed (pos & neg) elements and/or rows not allowed", 
#                              call. = FALSE)
#    if(is.numeric(j) && any(j<=0))
#        if(all(j<=0))
#            j <- c(1:ncol(x$data))[!1:ncol(x$data) %in% abs(j)] else
#                if(force=="omit.err.cases") 
#                     j <- j[j>0] else 
#                         stop("In pems[i,j]: mixed (pos & neg) elements and/or rows not allowed", 
#                              call. = FALSE)

    #make dummy data.frame for force == na.pad.output
    if(force=="na.pad.output"){
        dummy <- as.data.frame(matrix(nrow=length(i), ncol=length(j)), row.names = NULL, optional = TRUE)
        if(is.character(j)) names(dummy) <- j 
        if(is.numeric(j)) names(dummy) <- names(x$data)[j]
        row.names(dummy) <- i
    }

#try allows pems[does.not.exist,] 
#generates as data.frame of NAs
#so using another check

#note 
#might have fixed this later

#this need thinking about 
#could be rationalised


    index.i <- if(is.character(i))
                   !i %in% row.names(x$data) else !i %in% 1:nrow(x$data)
    check.i <- i[index.i]
    index.i <- c(1:length(index.i))[!index.i]

    index.j <- if(is.character(j))
                   !j %in% names(x$data) else !j %in% 1:ncol(x$data)
    check.j <- j[index.j]
    index.j <- c(1:length(index.j))[!index.j]

#old version -with no index.j and .j

##    check.i <- if(is.character(i))
##                   i[!i %in% row.names(x$data)] else i[!i %in% 1:nrow(x$data)]
##    check.j <- if(is.character(j))
##                   j[!j %in% names(x$data)] else j[!j %in% 1:ncol(x$data)]
  
    if(length(check.i)>0 || length(check.j)>0){

        #previously told user what was missing 
        #but that was messy if lots missing

        if("no.forcing" %in% force){
            temp <- paste(c("elements", "rows")[c(length(check.j)>0, length(check.i)>0)], sep="", collapse=" and ")
            stop(paste("In pems[i,j]: unknown ", temp, " called", sep="", collapse=""), 
                 call. = FALSE)
        }

        #these are know known i and j terms
        i <- i[!i %in% check.i]
        j <- j[!j %in% check.j]

        if("omit.err.cases" %in% force | "na.pad.output" %in% force){
            if(length(i)<1 || length(j)<1){
                temp <- paste(c("elements", "rows")[c(length(j)<1, length(i)<1)], sep="", collapse=" or ")
                if("omit.err.cases" %in% force)
                    stop(paste("In pems[i,j]: no known ", temp, " even after forcing", sep="", collapse=""), 
                         call. = FALSE)
             }
         }
        
    }

    #try to get x$data[i,j]
    ans <- try(x$data[i,j,..., drop=F], silent=TRUE)
    if(class(ans)[1]=="try-error" || (is.data.frame(ans) && nrow(ans)==0))
        if(force=="na.pad.output") 
            ans <- dummy else
                   stop("In pems[i,j], unexpected issue [please contact package admin.]", 
                        call. = FALSE)

#this is still in progress
#put everything that is there and then
#transfer the attributes?

    if("na.pad.output" %in% force){
         dummy[index.i,index.j] <- ans

    

##this works (I think...)
##but is well messy

##also suspect the following might die
##if units are present but does not contain
##an expected entry

        if(length(j)>0)
            for(jj in 1:length(j))
                attributes(dummy[,index.j[jj]]) <- attributes(ans[,jj])

        ans <- dummy
        dummy <- dummy[1,,drop = FALSE]
        dummy[index.j] <- x$units[j]
        names(ans) <- make.names(names(ans), unique = TRUE)
        names(dummy) <- names(ans)
        j <- names(ans)
        x$units <- dummy

#known issue
#this makes no unit cases NA
#elsewhere units "" if not set
#does this matter 
#which to do, if doing only one?
#if doing both, are both handled elsewhere?
    }

    #if simplify can be done return pems.element
    if(simplify & ncol(ans)==1){

        out <- ans[,1] 
        attr(out, "row.names") <- NULL
        attr(out, "name") <- names(ans)
        attr(out, "units") <- as.character(x$units[1,j])
        class(out) <- unique(c("pems.element", class(out)))
################
#think about makePEMSElement
################
        return(out)
    }
    
    #otherwise return rebuilt pems
    x$data <- ans
    x$units <- x$units[1,j]
    if("history" %in% names(x))
         x$history <- c(x$history, call2)
      
    class(x) <- old.class
    x

}




#########################
#########################
##[<-.pems
#########################
#########################

`[<-.pems` <- function(x, i, j, ..., force = FALSE, value){

    ########################
    #generic pems handling
    ########################

#overwrite not currently doing anything

#?option for pems[1,1] <- 1:10 to force insert as pems[1:10,1]
#even if pems[1:10,1] exists 
#this would be the overwrite???

#?option for pems[1:2] <- vector shorter than pems rows
#to write as element replicated as vector+NAs
#note: this is na.pad.target for elements then fill.target/insert for rows. 

    #x[1] element 1 as 1 col data.frame
    #x[1,] row 1
    #x[,1] element 1 
    #x[1,1] row, element 1,1

    #x[,] <- value insert value into  
    
    #output options
    #force=FALSE fail call if any requested elements/rows dimensions 
    #            do not fit exactly or is any missing 
    #force=TRUE return valid element/row requests, discard invalid 
    #            trim dimensions to trim smallest
    # 
#???    #simplify=TRUE as pems.element if possible else pems
#???    #simplify=FALSE as pems 

    #currently handles negs like data.frame (I think)...

    #negs are discarded if a mixed (neg and pos) vector supplied
    #and forced. Possibly a better way of handling this

#????    #overwrite = TRUE, if pems info/attributes in value 
    #                  write them over what is pems else use 
    #                  pems attributes if there
#????    #overwrite = FALSE, if pems info/attributes in pems 
    #                  write retain it. If nothing and 
    #                  something in value should be use it?
    #                  or does this make third case

#????    #attribute.source = "x", x only "x.value" x then value, etc, ...
    #                   "value", "value.x" 


    #######################
    #check what I was supplied
    #######################

    call1 <- sys.call()
    call2 <- match.call()
    check.i <- !missing(i)
    check.j <- !missing(j)
    check.op <- FALSE

    ###########################
    #force special handling
    ###########################

    #options omit.err.cases (default forcing)
    if(is.logical(force))
        force <- if(force) "omit.err.cases" else "no.forcing"
    if(!is.character(force))
         stop("In pems[i,j]<-: unknown force option", 
                              call. = FALSE)

# this would restrict forcing to known types

#    if(!all(force %in% c("no.forcing", "omit.err.cases", "na.pad.target", "na.pad.insert", "fill.target", "crop.insert")))
#         stop("In pems[i,j]<-: unknown force option", 
#                              call. = FALSE)

    ##############################
    #crack open pems
    ##############################
    
    old.class <- class(x)
    class(x) <- "not.pems"

    #########################
    #make it a standard pems[i,j]
    #########################

    #request pems[], pems[,], etc
    #might think about this 
    if(!check.i & !check.j) {
        i <- 1:nrow(x$data)
        j <- 1:ncol(x$data)
    }

    
    #request pems[i], etc
    #columns/elements, etc
    if(check.i && !check.j && length(as.character(call1))==length(as.character(call2))){
       j <- i
       i <- 1:nrow(x$data)
       check.op <- TRUE
    }

    #request pems[i,], etc
    #rows, etc
    if(check.i && !check.j && length(as.character(call1))!=length(as.character(call2)))
        j <- 1:ncol(x$data)

    #request pems[,j], etc
    #columns/elements
    if(!check.i && check.j)
        i <- 1:nrow(x$data)
    
    #at this point the request must be in form pems[i,j]

    ##############################
    #special case - negatives
    ##############################
 
    #current not going to accept them
    test <- character()
    if(is.numeric(j) && any(j<=0)) test <- c(test, "elements")
    if(is.numeric(i) && any(i<=0)) test <- c(test, "rows")
    if(length(test)>0)
        stop(paste("In pems[i,j]<-: neg ",
                   paste(test, sep=" and ", collapse=" and "),
                   " requested \n       [not currently allowed]",
                   sep=""), call. = FALSE)

    ################################
    #identify allowed cases
    ################################
    #T/F - is/isnt there 
    check.i <- if(is.character(i))
                   i %in% row.names(x$data) else i %in% 1:nrow(x$data)
    check.j <- if(is.character(j))
                   j %in% names(x$data) else j %in% 1:ncol(x$data)

    ###################
    #get value dimensions
    ###################

    #previous attempts to standardise value did not work
    check.value <- NULL
    if(is.vector(value) | is(value)[1]=="pems.element" | is.factor(value)){
        check.value <- "vector"
        value.dim <- c(length(value),1)
    }
    if(is.data.frame(value)){
        check.value <- "data.frame"
        value.dim <- dim(value)
    }
    if(is(value)[1]=="pems"){
        check.value <- "pems"
        class(value) <- "not.pems"
        pems.units <- value$units
        value <- value$data
        value.dim <- dim(value)
    }
    if(is.null(check.value)){
        stop("In pems[i,j]<-: can't insert insert of that class!", call.=FALSE)
    }

    

#    if(is.vector(value) | is(value)[1]=="pems.element"){
#        temp <- attributes(value)
#        value <- as.data.frame(value, drop=F, stringsAsFactors=FALSE)
#        if(!is.null(temp)) attributes(value[,1]) <- temp 
#    }
#    if(is.data.frame(value)){
#replace with makePEMS???
#need to make sure 1x1 data.frame 
#does not revert to vector
#        value <- list(data=as.data.frame(value, drop=F, stringsAsFactors =FALSE),
#                      units=data.frame())
#        class(value)<- "not.pems"
#    } 
#    if(class(value)[1]=="pems"){
#        class(value)<- "not.pems"
#    }
#    if(class(value)[1]!="not.pems"){
#        stop("In pems[i,j]<-: can't insert insert of that class!", call.=FALSE)
#    }

#    #here value must have structure value$data
    
    #######################
    #get dim of value 
    #######################
#    value.dim <- dim(value$data)

    #value.dim = value[rows, elements]

#################
#fixed for update
#################

    #######################
    #force by omit.err.cases
    #######################

    if("omit.err.cases" %in% force){
        i <- i[check.i]
        j <- j[check.j]
        check.i <- check.i[check.i]
        check.j <- check.i[check.j]
    }


#################
#fixed in update
#################

    #########################
    #force by na.pad.target (1)
    ######################### 

    #second part to this after crop.value

    if("na.pad.target" %in% force){
        if(any(!check.j)){
            if(is.character(j)){
                temp.j <- j[!check.j]
                j <- j[check.j]
                temp <- names(x$data)
                x$data[,temp.j] <- NA
                temp.j <- names(x$data)[!names(x$data) %in% temp]
#this could still fall over
#if data and unit names mismatched before
                x$units[,temp.j] <- NA
                j <- c(j, temp.j)
                check.j <- rep(TRUE, length(j))
            } else {
                #assuming it is numeric
                temp.j <- (ncol(x$data)+1) : max(j[!check.j], na.rm=TRUE)
                temp <- names(x$data)
                x$data[,temp.j] <- NA
                temp.j <- names(x$data)[!names(x$data) %in% temp]
#as above but remember 
#,10 would create 10 
#plus missing before it
                x$units[,temp.j] <- NA
                check.j <- rep(TRUE, length(j))
            }
        }
        #like above but for row.name, nrow, i, etc...
        #no unit update
        if(any(!check.i)){
            if(is.character(i)){
                temp.i <- i[!check.i]
                i <- i[check.i]
                temp <- row.names(x$data)
                x$data[temp.i,] <- NA
                temp.i <- row.names(x$data)[!row.names(x$data) %in% temp]
                i <- c(i, temp.i)
                check.i <- rep(TRUE, length(i))
            } else {
                #assuming it is numeric
                temp.i <- (nrow(x$data)+1) : max(i[!check.i], na.rm=TRUE)
                temp <- row.names(x$data)
                x$data[temp.i,] <- NA
                check.i <- rep(TRUE, length(i))
            }
        }        
    }

###################
#working for update
###################

    #########################
    #force by crop.insert
    ######################### 

    if("crop.insert" %in% force){
       if(value.dim[1] > length(i)){
           if(check.value=="vector")
               value <- rep(value, length.out=value.dim[1])
           if(check.value=="data.frame" | check.value=="pems")
               value <- value[1:length(i),]
           value.dim[1]<- length(i)
       } 
       if(value.dim[2] > length(j)){
           #should not happen?? with vector
           if(check.value=="data.frame")
               value <- value[,1:length(j)]
           if(check.value=="pems"){
               value$data <- value$data[,1:length(i)]
               pems.units <- pems.units[,1:length(i)]
           }
       }
    }

#######################
#fixed in update
#######################

    #########################
    #force by na.pad.target part 2
    #########################

    #this must be after crop or we pad for something in value
    #we then remove...

    #just do these as numerics for now
    #could look into assigning names from value if there later???

    if("na.pad.target" %in% force){
       if(length(j) < value.dim[2]){ 
           #assuming it is numeric
           temp.j <- (ncol(x$data)+1) : (ncol(x$data)+(value.dim[2]-length(j)))
           temp <- names(x$data)
           x$data[,temp.j] <- NA
           temp.j <- names(x$data)[!names(x$data) %in% temp]
#as na.pad.target 1 
           x$units[,temp.j] <- NA
           j <- c(j, temp.j)
           check.j <- rep(TRUE, length(j))
       }
       if(length(i) < value.dim[1]){ 
           temp.i <- (nrow(x$data)+1) : (nrow(x$data)+(value.dim[1]-length(i)))
           temp <- row.names(x$data)
           x$data[temp.i,] <- NA
           temp.i <- row.names(x$data)[!row.names(x$data) %in% temp]
           i <- c(i, temp.i)
           check.i <- rep(TRUE, length(i)) 
       }
    }

#####################
#fixed for update
#####################

    #############################
    #na.pad.insert
    #############################
    if("na.pad.insert" %in% force){
        if(length(j) > value.dim[2]){
            #if value is a vector don't need to do anything
            #except tell it that it is OK
            if(check.value=="vector"){
                value.dim[2] <- length(j) 
            } else {
                temp.j <- (ncol(value)+1) : (ncol(value)+(length(j)-value.dim[2]))
                value[,temp.j] <- NA
                value.dim <- dim(value)
            }
        }
        if(length(i) > value.dim[1]){
            if(check.value=="vector"){
                value <- c(value, rep(NA, length.out=length(i) - value.dim[1]))
                value.dim[1] <- length(i)
            } else {
                temp.i <- (nrow(value)+1) : (nrow(value)+(length(i)-value.dim[1]))
                value[temp.i,] <- NA            
                value.dim <- dim(value)
            }
        }
    }


#####################
#fixed for update
#####################

    #############################
    #fill.insert
    #############################
    if("fill.insert" %in% force){
        if(length(i) > value.dim[1]){
            if(check.value=="vector"){
                value <- c(value, rep(value, length.out=length(i) - value.dim[1]))
                value.dim[1] <- length(i)
            } else {
                temp.i <- (nrow(value)+1) : (nrow(value)+(length(i)-value.dim[1]))
                value[temp.i,] <- value[rep(1:nrow(value), length.out=length(temp.i)),]            
                value.dim <- dim(value)
            }
        }
        if(length(j) > value.dim[2]){
            #if value is a vector don't need to do anything
            #except tell it that it is OK
            if(check.value=="vector"){
                value.dim[2] <- length(j) 
            } else {
                temp.j <- (ncol(value)+1) : (ncol(value)+(length(j)-value.dim[2]))
                value[,temp.j] <- value[,rep(1:ncol(value), length.out=length(temp.j))]
                value.dim <- dim(value)
            }
        }
    }

    ###############################
    #check dimensions
    ###############################

#could move the first of these message to just after 
#omit.err.cases
#might fall over otherwise

    #error out if it does not fit.
    fault.message <- "In pems[i,j]<-value: questionable request"
    if(length(i)<1 | length(j)<1){
        temp <- paste(c("elements", "rows")[c(length(j)<1, length(i)<1)], sep="", collapse=" and ")
        fault.message <- paste(fault.message, "\n       no valid ", temp, " set in pems[i,j]", sep="", collapse="")
        fault.message <- paste(fault.message, "\n       [check force setting if insertion required]", sep="", collapse="")
        stop(fault.message, call.=FALSE)
    }
    if(any(!check.i) | any(!check.j)){
        temp <- paste(c("elements", "rows")[c(any(!check.j), any(!check.i))], sep="", collapse=" and ")
        fault.message <- paste(fault.message, "\n       unknown ", temp, " set in pems[i,j]", sep="", collapse="")
        fault.message <- paste(fault.message, "\n       [check force setting if insertion required]", sep="", collapse="")
        stop(fault.message, call.=FALSE)
    }
    if(length(i) != value.dim[1] | length(j) != value.dim[2]){
        temp <- paste(c("elements", "rows")[c(length(j)!=value.dim[2], length(i)!=value.dim[1])], sep="", collapse=" and ")
        fault.message <- paste(fault.message, "\n       ", temp, " pems[i,j] and insert[i,j] dimension mismatch", sep="", collapse="")
        fault.message <- paste(fault.message, "\n       [check force setting if insertion required]", sep="", collapse="")
        stop(fault.message, call.=FALSE)
    }

    #################################
    #try insert
    #################################

    test <- if(check.op)
                try(x$data[j] <- value, silent=T) else 
                try(x$data[i,j] <- value, silent=T)
    if(class(test)[1]=="try-error")
        stop("In pems[i,j]<-value: bad insertion \n       significant pems[i,j]/value class mismatch",
              call.=FALSE)



#error messaging could be more informative?
#could like this into addition force/overwrite options

    ##################################
    #update units
    ##################################
    
    #currently
    #only copying attributes if units not there
    #options could be possible, overwrite or force options
    #could also do a method for data.frames or pems without units
    #would read attributes of each column

    if(check.value=="vector"){
        if("units" %in% names(attributes(value))){
            for(jj in j){
######################
#update units of input 
#replaces units of pems
#######################
##               if(check.op || is.null(x$units[1,jj]) || is.na(x$units[1,jj]) || x$units[1,jj]=="")
                  x$units[1,jj] <- attributes(value)$units
##                  attributes(x$data[,jj]) <- attributes(value)
               }
        }
    }    
    if(check.value=="pems"){
#might need an error catcher for 
#not pems.units
         for(jj in 1:length(j))
######################
#update units of input 
#replaces units of pems
#######################
##             if(is.null(x$units[1,j[jj]]) || is.na(x$units[1,j[jj]]) || x$units[1,j[jj]]==""){
                 x$units[1,j[jj]] <- pems.units[1,jj]
##                 attributes(x$data[,j[jj]]) <- attributes(value[,jj])
##             }
    }    

    #####################
    #history update
    #####################

    #might want a silence history logging option?

    if ("history" %in% names(x)) 
        x$history <- c(x$history, call2)

################################
#class conflicts not currently handled
#next jobs....
#think about simplify
#does it have meaning here?
################################

    ####################################
    # send back data
    ####################################

    class(x) <- old.class
    return(x)

}











##########################
#########################
##$.pems
#########################
#########################

`$.pems` <- function(x, name, ...){


#####################
#old version
#####################
#    class(x) <- "not.pems"
#    ans <- try(x$data[, i], silent = TRUE)
#    if(class(ans)[1] == "try-error"){
#        warning("Element '", i, "' not found in pems", call. = FALSE)
#        return(NULL)
#    }
#    if (!is.null(ans)) 
#        attr(ans, "name") <- i
#    if (!is.null(ans) && !is.null(units)) 
#        if (is.null(attributes(ans)$units)) 
#            attr(ans, "units") <- x$units[1,i]
#        class(ans) <- "pems.element"
#    ans
#think about x[,i, simplify=TRUE]
#might not be need because

#######################
#another old version
#######################
#    ans <- try(x[, i], silent = TRUE)
#    if(class(ans)[1] == "try-error"){
#        warning("Element '", i, "' not found in pems", call. = FALSE)
#        return(NULL)
#    }
#    if (!is.null(ans)) 
#        attr(ans, "name") <- i
##########################
#this looks a bit screwy
#########################
#    if (!is.null(ans) && !is.null(units)) 
#        if (is.null(attributes(ans)$units)) 
#            attr(ans, "units") <- x$units[1,i]
#    ans

########################
#another old version
########################
#    ans <- try(x[, name, simplify = TRUE, force = FALSE], silent = TRUE)
#make silent like data.frame$does.not.exist
#    if(class(ans)[1] == "try-error"){
#        warning("Element '", i, "' not found in pems", call. = FALSE)
#        return(NULL)
#    } else return(ans)

    ans <- try(x[, name, simplify = TRUE, force = FALSE], silent = TRUE)
    if(class(ans)[1] == "try-error") NULL else ans

}



#########################
#########################
##$<-.pems
#########################
#########################

`$<-.pems`<- function(x, name, ..., value){

      x[name, force=c("na.pad.insert", "na.pad.target")] <- value
      x

}





#########################
#########################
##with.pems
#########################
#########################

#version 0.1.0 kr 2015-08-02

#test run parsync/carb data analysis

#note currently discards units
#not sure there is a way around this...

with.pems <- function(data, expr, ...) {

   eval(substitute(expr), pemsData(data), enclos = parent.frame())

}





#########################
#########################
##subset.pems
#########################
#########################

#version 0.1.0 kr 2015-09-07

#test run parsync/ucr data analysis

subset.pems <- function(x,...){

    x[["data"]] <- subset(x[["data"]], ...)
    x

}






##########################
##########################
##print.pems
##########################
##########################

#kr 06/06/2013 v 0.3.0

#what it does
##########################
#handles pems console appearance 
#etc
#

#to do
##########################
#remove temp?


print.pems <- function (x, verbose = FALSE, n=6, ...) {

    temp <- x
    class(temp) <- "not.pems"

    #show data frame

###########################
#think about this and head
###########################

    if(verbose){
        print.data.frame(temp$data)
    } else {


        if(is.null(nrow(temp$data)) || nrow(temp$data)<=n) 
           print.data.frame(temp$data) else {
           print.data.frame(temp$data[1:n, , drop=FALSE])
           message("...")
        } 
    }

    #pems report line 1

    reply <- names(temp$data)
    if (is.null(reply)) 
        message("\npems object: no named data [suspect]")
    else message("\npems object: ", ncol(temp$data), " data series (each ", 
        nrow(temp$data), " cases)")

    #pems report line 2 structure

    reply <- names(temp)[names(temp) %in% c("units", "constants", 
        "history")]
    if (length(reply) < 1) 
        message("\twith no supporting structure [suspect]")
    else message("\twith supporting structures: ", paste(reply, 
        collapse = ", ", sep = ""))

    #pems report line 3 extra tags

    reply <- names(temp)[!names(temp) %in% c("data", "units", 
        "constants", "history", "dem")]
    if (length(reply) > 0) 
        message("\t[and unique tags: ", paste(reply, collapse = ", ", 
            sep = ""), "]\n")

    #output
    invisible(x)

}



##########################
##########################
##plot.pems
##########################
##########################

#kr 07/12/2011 v 0.2.0

#what it does
##########################
#generates simple plot
#does not keep units
#

#to do
##############################
#remove temp


##' @S3method plot pems
plot.pems <- function(x, id = NULL, ignore = "time.stamp", n = 3, ...) {

   temp <- x
   class(temp) <- "no.class"

   reply <- temp$data

   if(is.null(reply)){
      message("\npems object [suspect]")
      return(invisible(NULL))
   }

   if(is.null(id)){
       id <- 1:ncol(reply)
       if(length(ignore)>0)
           id <- id[!names(reply) %in% ignore]
       if(n>0 && length(id)>n)
           id <- id[1:n]
   }

   reply <- reply[id]      

   plot(reply, ...)

}



##########################
##########################
##names.pems
##########################
##########################

#kr 07/12/2011 v 0.2.0

#what it does
##########################
#returns data series names from pems

#to do
##########################
#names<- handling
#



##' @S3method print pems
names.pems <- function(x, ...) {

    class(x) <- "no.class"
    x <- x$data

    if(is.null(x)){
       message("\npems object [suspect]")
       return(invisible(NULL))
    }

    names(x)
}


## @S3method names<-.pems

`names<-.pems` <- function(x, ..., value) {

    #variation on units<-
    #very crude handling of names$ and names[]

    call2 <- match.call()
    class(x)[1] <- "no.class"

#    if(***************){
#       warning("In units(pems): [suspect units]", call.=FALSE)
#       units <- data.frame(matrix(NA, nrow = 1, ncol = ncol(x$data)))
#       names(units) <- names(x$data)
#    }

#might also think about adding units if not there
#see above

#error check on state of value?

    names(x$data) <- value
    names(x$units) <- value

    if ("history" %in% names(x)) 
        x$history <- c(x$history, call2)
    
    class(x)[1] <- "pems"
    return(x)

}







##########################
##########################
##summary.pems
##########################
##########################

#kr 07/12/2011 v 0.2.0

#what it does
##########################
#generates summary reports 
#

#to do
##########################
#make dedicated summary 
#and option for as current as alternative
#


#comments
##########################
#to tidy


##' @S3method print pems
summary.pems <- function(object, ...) {

   class(object) <- "no.class"

   object <- object$data

   if(is.null(object)){
      message("\npems object [suspect]")
      return(invisible(NULL))
   }

   return(summary(object))
                       
}





##########################
##########################
##units.pems
##########################
##########################

#kr 07/12/2011 v 0.2.0

#what it does
##########################
#extracts units from pems 
#

#to do
##########################
#make dedicated summary 
#and option for as current as alternative
#


#comments
##########################
#to tidy


##' @S3method units.pems
units.pems <- function(x) {

    class(x) <- "no.class"
    x <- x$units

    if(is.null(x)){
       warning("In units(pems): pems unitless [suspect]", call.=FALSE)
       return(invisible(NULL))
    }

    return(x)

}


## @S3method unit<-.pems

`units<-.pems` <- function(x, value) {

    call2 <- match.call()
    class(x)[1] <- "no.class"
    units <- x$units

    if(is.null(units)){
       warning("In units(pems): [suspect units]", call.=FALSE)
       units <- data.frame(matrix(NA, nrow = 1, ncol = ncol(x$data)))
       names(units) <- names(x$data)
    }

#should this be a try?
#should this foreshorten/expand?
#should be check bad moves

#######################

#way units methods works this is only every a data.frame and what you have to send it have to work
#so thinking what follows is a waste of time
#could just be x$units <- value

    if(is.data.frame(value)){

       units[1, 1:ncol(value)] <- value[1,]
       if(!is.null(names(units)))
          names(units)[1:ncol(value)] <- names(value)
    } else {

       units[1,1:length(value)] <- as.character(value)
    }

########################


    if ("history" %in% names(x)) 
        x$history <- c(x$history, call2)
    
    x$units <- units
    class(x)[1] <- "pems"
    return(x)

}




##########################
##########################
##head.pems
##tail.pems
##########################
##########################

#kr 31/04/2014 v 0.2.4

#if number of columns less than n?
#need to catch this


head.pems <- function(x, n=6, ...){
    out <- x[1:n,,force=T,simplify=F]
    print(out, verbose=TRUE)
    invisible(out)
}


tail.pems <- function(x, n=6, ...) {
    out <- dim(as.data.frame(x))[1]
    out <- (out-n+1):(out)
    out <- out[out>0]
    out <- x[out,,force=TRUE,simplify=F]
    print(out, verbose=TRUE)
    invisible(out)
}







############################################################
############################################################








###############################
#working
###############################



###########################
#i do not think you can do anymore than
#units(pems) <- value
##########################
#error out if it does not work
#not even sure we can do this
##########################



