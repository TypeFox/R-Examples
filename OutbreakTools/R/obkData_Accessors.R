
###############################
#### ACCESSORS FOR OBKDATA ####
###############################

###############
## get.locus ##
###############
setMethod("get.locus", "obkData", function(x, ...){
    if(is.null(x@dna)) return(NULL)
    return(get.locus(x@dna, ...))
})



################
## get.nlocus ##
################
setMethod("get.nlocus", "obkData", function(x, ...){
    if(is.null(x@dna)) return(0)
    return(get.nlocus(x@dna, ...))
})




###################
## get.sequences ## (get sequence ID)
###################
setMethod("get.sequences","obkData", function(x, ...){
    return(get.sequences(x@dna, ...))
})




####################
## get.nsequences ##
####################
setMethod("get.nsequences", "obkData", function(x, ...){
    if(is.null(x@dna)) return(0)
    return(get.nsequences(x@dna, ...))
})



#############
## get.dna ##
#############
setMethod("get.dna", "obkData", function(x, locus=NULL, id=NULL, ...){
    ## checks and escapes ##
    if(is.null(x@dna)) return(NULL)
    return(get.dna(x@dna, locus=locus, id=id, ...))
})



#####################
## get.individuals ##
#####################
setMethod("get.individuals", "obkData", function(x, data=c("all", "individuals", "records", "contacts", "dna"), ...){
    data <- match.arg(data)

    ## list individuals in @individuals
    if(data=="individuals"){
        if(is.null(x@individuals)) return(NULL)
        return(row.names(x@individuals))
    }

    ## list individuals in @records
    if(data=="records"){
        if(is.null(x@records)) return(NULL)
        out <- unlist(lapply(x@records, function(e) e$individualID), use.names=FALSE)
        return(unique(out))
    }

    ## list individuals in @contacts
    if(data=="contacts"){
        if(is.null(x@contacts)) return(NULL)
        return(get.individuals(x@contacts))
    }

    ## list individuals in @contacts
    if(data=="dna"){
        if(is.null(x@dna)) return(NULL)
        return(get.individuals(x@dna))
    }


    ## list all individuals in the object (in @individuals, @records and @contacts)
    if(data=="all"){
        out <- unlist(lapply(c("individuals", "records", "contacts", "dna"), function(e) get.individuals(x, e)), use.names=FALSE)
        return(unique(out))
    }
})



######################
## get.nindividuals ##
######################
setMethod("get.nindividuals", "obkData", function(x, data=c("all", "individuals", "records", "contacts", "dna"), ...){
    data <- match.arg(data)

    return(length(get.individuals(x, data=data)))
})




#################
## get.records ##
#################
setMethod("get.records", "obkData", function(x, ...){
  if(is.null(x@records)) return(NULL)
  return(names(x@records))
})



##################
## get.nrecords ##
##################
setMethod("get.nrecords", "obkData", function(x, ...){
  if(is.null(x@records)) return(0)
  return(length(get.records(x)))
})


#################
## get.context ##
#################
setMethod("get.context", "obkData", function(x, ...){
  if(is.null(x@context)) return(NULL)
  return(names(x@context))
})



##################
## get.nrecords ##
##################
setMethod("get.ncontext", "obkData", function(x, ...){
  if(is.null(x@context)) return(0)
  return(length(get.context(x)))
})


###############
## get.dates ##
###############
setMethod("get.dates", "obkData", function(x, data=c("all", "individuals", "records", "dna", "context", "contacts"),...){

  data <- match.arg(data)

    ## list dates in @individuals
    if(data=="individuals"){
        if(is.null(x@individuals)) return(NULL)
        return(x@individuals$date)
    }

    ## list dates in @records
    if(data=="records"){
        if(is.null(x@records)) return(NULL)
        out <- unlist(lapply(x@records, function(e) as.character(e$date)))
        return(as.Date(out))
    }

    ## list dates in @context
    if(data=="context"){
      if(is.null(x@context)) return(NULL)
      out <- unlist(lapply(x@context, function(e) as.character(e$date)))
      return(as.Date(out))
    }

    ## list individuals in @dna
    if(data=="dna"){
        if(is.null(x@dna)) return(NULL)
        return(get.dates(x@dna))
    }

  ## list individuals in @contacts
    if(data=="contacts"){
        if(is.null(x@contacts)) return(NULL)
        return(get.dates(x@contacts))
    }


    ## list all individuals in the object (in @individuals, @records and @contacts)
    if(data=="all"){
        out <- unlist(lapply(c("individuals", "records", "dna"), function(e) as.character(get.dates(x, e))))
        return(as.Date(out))
    }
})



################
## get.ndates ##
################
setMethod("get.ndates", "obkData", function(x, data=c("all", "individuals", "records", "dna","context", "contacts"),...){
    return(length(get.dates(x, data=data)))
})



###############
## get.ntrees ##
###############
setMethod("get.ntrees", "obkData", function(x, ...){
    if(is.null(x@trees)) return(0L)
    return(length(x@trees))
})


###############
## get.trees ##
###############
setMethod("get.trees", "obkData", function(x, ...){
    return(x@trees)
})


##################
## get.contacts ##
##################
setMethod("get.contacts", "obkData", function(x, from=NULL, to=NULL, ...){
    if(is.null(x@contacts)) return(NULL)
    return(get.contacts(x@contacts, from=from, to=to, ...))
})



###################
## get.ncontacts ##
###################
setMethod("get.ncontacts", "obkData", function(x, from=NULL, to=NULL, ...){
    if(is.null(x@contacts)) return(0)
    return(get.ncontacts(x@contacts, from=from, to=to, ...))
})



##############
## get.data ##
##############
##
## Universal accessor:
## tries to find any type of data within the obkData object
##
setMethod("get.data", "obkData", function(x, data, where=NULL, drop=TRUE, showSource=FALSE, ...){
    ## disable bloody stringsAsFactor ##
    o.saf <- options("stringsAsFactors")
    on.exit(options(o.saf))
    options(stringsAsFactors=FALSE)

    data <- as.character(data)

    result <- data.frame()

    ## LOOK FOR SLOT NAMES ##
    if(data[1] %in% slotNames(x)) return(slot(x, data))

    ## HANDLE 'WHERE' ##
    if(!is.null(where)){
        where <- match.arg(as.character(where), c("individuals", "records", "context", "dna"))

        ## look in @individuals ##
        if(where=="individuals"){
            if(is.null(x@individuals)) { # return NULL if empty
                warning("x@individuals is NULL")
                return(NULL)
            }
            if(any(data %in% names(x@individuals))){
                temp<-x@individuals[,data,drop=FALSE] # get data
                temp<-cbind(temp,rownames(x@individuals)) # add individualID
                temp<-cbind(temp,rep("individuals",dim(temp)[1])) # add source
                result<-temp
                names(result)<-c(data,"individualID", "source")
            } else {
                warning(paste("data '", data, "' was not found in @individuals",sep=""))
                return(NULL)
            }
        } # end where==individuals


        ## look in @records ##
        if(where=="records"){
            if(is.null(x@records)) { # return NULL if empty
                warning("x@records is NULL")
                return(NULL)
            }
            ## look in @records
            if(any(data %in% names(x@records))){
                return(x@records[[data]])
            }
            ## look within slots in @records
            found <- FALSE
            for(i in 1:length(x@records)){
                if(any(data %in% names(x@records[[i]]))){
                    found <- TRUE
                    temp<-x@records[[i]][,c(data,"individualID","date")]
                    temp<-cbind(temp,rep(names(x@records)[i],dim(temp)[1]))
                    colnames(temp)<-c(data,"individualID","date","source")
                    result<-rbind(result,temp)
                }
            }
            if(!found){
                warning(paste("data '", data, "' was not found in @records", sep=""))
                return(NULL)
            }
        } # end where==records


        ## look in @context ##
        if(where=="context"){
            if(is.null(x@context)) { # return NULL if empty
                warning("x@context is NULL")
                return(NULL)
            }
            ## look in @context
            if(any(data %in% names(x@context))){
                return(x@context[[data]])
            }
            ## look within elements of @context
            found <- FALSE
            for(i in 1:length(x@context)){
                if(any(data %in% names(x@context[[i]]))){
                    found <- TRUE
                    temp<-x@context[[i]][,c(data,"date")]
                    temp<-cbind(temp,rep(names(x@context)[i],dim(temp)[1]))
                    colnames(temp)<-c(data,"date","source")

                    result<-rbind(result,temp)
                }
            }
            if(!found){
                warning(paste("data '", data, "' was not found in @context",sep=""))
                return(NULL)
            }
        } # end where==context

        ## look within elements of @dna ##
        if(where=="dna"){
            if(is.null(x@dna)){ # return NULL if empty
                warning("x@dna is NULL")
                return(NULL)
            }
            if(any(data %in% names(x@dna@meta))){
                temp <- x@dna@meta[, c(data, "individualID", "date"), drop=FALSE]
                temp <- cbind(temp,rep("dna",dim(temp)[1]))
                result <- temp
                names(result) <- c(data,"individualID", "date", "source")
            } else {
                warning(paste("data '", data, "' was not found in @dna",sep=""))
                return(NULL)
            }
        } # end where==dna


    } # end if 'where' provided
    else{
        ## ELSE, LOOK EVERYWHERE ##
        result <- lapply(c("individuals", "records", "context", "dna"), function(e)
                         suppressWarnings(get.data(x, data=data, where=e, drop=drop, showSource=TRUE)))

        result <- Reduce("rbind", result)
        if(!is.null(result) && nrow(result)>0) rownames(result) <- NULL # rownames are meaningless

    } # end search everywhere


    ## RETURN REQUESTED RESULT ##
    if(length(result)>0){
        if(showSource)
            return(result)
        else
            if(any(data %in% names(result))){
                data <- data[data %in% names(result)]
                result <- result[,data,drop=drop]
            }
            return(result)
    }
    else{
        ## DEFAULT IF WE DON'T KNOW WHAT TO RETURN ##
        warning(paste("data '", data, "' was not found in the object", sep=""))
        return(NULL)
    }
}) # end get.data





