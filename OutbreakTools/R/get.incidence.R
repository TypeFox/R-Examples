

###################
## get.incidence ##
###################
setGeneric("get.incidence", function(x, ...) standardGeneric("get.incidence"))



#################
## Date method ##
#################
setMethod("get.incidence", "Date", function(x, from=NULL, to=NULL,
                                            interval=1, add.zero=TRUE, ...){
    ## GET DATES OF THE OUTPUT ##
    ## get first/last dates ##
    first.date <- min(x, na.rm=TRUE)
    last.date <- max(x, na.rm=TRUE)

    ## handle from ##
    if(is.null(from)){
        from <- first.date
    }
    if(is.numeric(from)) from <- first.date+from
    if(!inherits(from, "Date")) from <- as.Date(from)
    from <- min(from, first.date)

    ## handle to ##
    if(is.null(to)){
        to <- last.date
    }
    if(is.numeric(to)) to <- first.date+to
    if(!inherits(to, "Date")) to <- as.Date(to)
    to <- max(to, last.date)

    ## generate output dates ##
    interval <- round(interval)
    out.dates <- seq(from, to, by=interval) # output dates


    ## COMPUTE INCIDENCE ##
    ## incid is computed for time intervals
    ## incid on interval [d1,d2[ is named after d1
    breaks <- c(as.integer(out.dates), as.integer(to)+interval)
    incid <- table(cut(as.integer(x), breaks=breaks, right=FALSE))

    out <- data.frame(date=out.dates, incidence=as.integer(incid))

    ## add zero at the end if needed
    if(add.zero && incid[length(incid)] > 1e-14){
        out <- as.list(out)
        out$date <- c(out$date, to+interval)
        out$incidence <- c(out$incidence, 0L)
        out <- as.data.frame(out)
    }

    return(out)
}) # end Date method





#########################
## obkSequences method ##
#########################
setMethod("get.incidence", "obkSequences", function(x, from=NULL, to=NULL,
                                                    interval=1, add.zero=TRUE, ...){
    if(is.null(x) || get.nsequences(x)<1) return(NULL)

    out <- get.incidence(x@meta$date, from=from, to=to,
                         interval=interval, add.zero=add.zero, ...)
    return(out)

}) # end obkSequences method






########################
## obkContacts method ##
########################
setMethod("get.incidence", "obkContacts", function(x, from=NULL, to=NULL,
                                                   interval=1, add.zero=TRUE, ...){
    if(is.null(x) || get.ncontacts(x)<1 || !is.networkDynamic(x@contacts)) return(NULL)

    ## CHECK THAT THIS IS A DYNAMIC CONTACT NETWORK ##
    out <- get.incidence(as.data.frame(x)$onset, from=from, to=to,
                         interval=interval, add.zero=add.zero, ...)
    return(out)

}) # end obkContacts method




####################
## obkData method ##
####################
##
## based on 'dates' associated to a given field
## 'values' are optional and can be used to subset the retained 'dates'
## (e.g. define what a positive case is)
setMethod("get.incidence", "obkData", function(x, data, where=NULL, val.min=NULL, val.max=NULL, val.kept=NULL, regexp=NULL,
                                               from=NULL, to=NULL, interval=1, add.zero=TRUE, ...){
    ## HANDLE ARGUMENTS ##
    if(is.null(val.min)) val.min <- -Inf
    if(is.null(val.max)) val.max <- Inf


    ## GET DATA ##
    df <- get.data(x, data=data, where=where, showSource=TRUE)
    if(is.null(df)) stop(paste("Data",data,"cannot be found in this obkData object"))

    ## call specific procedures if applicable ##
    if(inherits(df, c("obkSequences", "obkContacts"))) {
        return(get.incidence(df, from=from, to=to,
                             interval=interval, add.zero=add.zero))
    }


    ## OTHERWISE: DATA ASSUMED TAKEN FROM RECORDS ##
    ## if data=='records', keep the first data.frame of the list ##
    if(is.list(df) && !is.data.frame(df) && is.data.frame(df[[1]])) df <- df[[1]]

    ## get dates ##
    if(!"date" %in% names(df)) stop("no date in the data")
    dates <- df$date

    ## get optional values associated to the dates ##
    ## keep 'data' if it is there
    if(data %in% names(df)){
        values <- df[[data]]
    } else { ## else keep first optional field
        temp <- !names(df) %in% c("individualID","date") # fields being not "individualID" or "date"
        if(any(temp)) {
            values <- df[,min(which(temp))]
        } else {
            values <- NULL
        }
    }


    ## EXTRACT RELEVANT DATES ##
    if(!is.null(values)){
        toKeep <- rep(TRUE, length(values))

        ## if 'values' is numeric ##
        if(is.numeric(values)){
            toKeep <- toKeep & (values>=val.min & values<=val.max)
        }

        ## if val.kept is provided ##
        if(!is.null(val.kept)) {
            toKeep <- toKeep & (values %in% val.kept)
        }

        ## if regexp is provided ##
        if(!is.null(regexp)) {
            temp <- rep(FALSE, length(values))
            temp[grep(regexp, values, ...)] <- TRUE
            toKeep <- toKeep & temp
        }

        dates <- dates[toKeep]
    }

    ## CALL THE DATE PROCEDURE ##
    if(length(dates)==0) return(NULL)
    out <- get.incidence(dates, from=from, to=to,
                         interval=interval, add.zero=add.zero)

    ## RETURN OUTPUT ##
    return(out)
}) # end obkData method

