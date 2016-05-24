


## hack to remove the NOTE in R CMD check about:
## plotIndividualTimeline: no visible binding for global variable ‘yITL’
## plotIndividualTimeline: no visible binding for global variable ‘value’
## plotIndividualTimeline: no visible binding for global variable ‘type’
if(getRversion() >= "2.15.1")  utils::globalVariables(c("yITL","date","type", "individualID"))



## Function to plot a timeline for individuals involved in an outbreak in one plot

.meltDateProof <- function(data,id.vars,measure.vars,variable.name){
    ## 'melt' in the reshape package looses the Date format, here is a hack around
    if(class(data[,measure.vars[1]])=="Date"){
        for(e in measure.vars){
            ## convert all columns (which should be dates) to chars
            data[,e] <- as.character(data[,e])
        }
        df <- melt(data,id.vars=id.vars,measure.vars=measure.vars,variable.name=variable.name)
        ## and convert back to Date
        df$value <- as.Date(df$value)
    }
    else
        df <- melt(data,id.vars=id.vars,measure.vars=measure.vars,variable.name=variable.name)
    return(df)
}


##########################
## plotIndividualTimeline
##########################
plotIndividualTimeline <- function(x, what="", selection=NULL, ordering=NULL, orderBy=NULL, colorBy=NULL,
                                   periods=NULL, plotNames=length(selection)<50, ...){
    ## plot selection of the individuals in data as a line, ordered by ordering
    ## color by colorby
    ## make lines for periods, an Nx2 matrix of column names

    if(is.null(selection)) selection <- 1:get.nindividuals(x, "individuals")
    if(length(selection)>length(get.individuals(x))){
        warning("selection is longer than the number of individuals. Selecting all.")
        selection <- 1:get.nindividuals(x, "individuals")
    }
    if(is.null(ordering)) ordering <- 1:length(selection)

    ## GET THE DATA.FRAME FOR THE TIME SERIES ##
    ## get subset of the plotted data ##
    x <- subset(x, individuals=selection)

    ## define ordering ##
    if(!is.null(orderBy)){
        ## alternatively, order by this character
        ordering <- order(unlist(get.data(x, data=orderBy, where="individuals")))
        ## ordering <- sapply(1:length(selection), function(i) which(order(get.data(x,orderBy)[selection])==i))
    }

    ## get time series data.frame ##
    df.ts <- make.individual.attributes(x)[ordering, , drop=FALSE]
    individualID <- df.ts$individualID

    ## isolates dates matching with requested data
    toKeep <- unlist(lapply(what, grep, names(df.ts)))
    if(length(toKeep)<1) stop(paste(what, "was not found in the data"))
    df.ts <- df.ts[, toKeep,drop=FALSE]

    ## keep only dates
    toKeep <- sapply(df.ts, inherits, "Date")
    if(length(toKeep)<1) stop("No date information to be used for the timeline plot")
    df.ts <- df.ts[, toKeep,drop=FALSE]

    ## convert dates to characters and back
    isDate <- sapply(df.ts, inherits, "Date")
    for(i in 1:ncol(df.ts)) if(inherits(df.ts[[i]], "Date")) df.ts[[i]] <- as.character(df.ts[[i]])

    ## add individualID and melt
    df.ts$individualID <- individualID
    df.ts <- melt(df.ts, id.var=c("individualID"))

    ## remove NAs, restore Date class
    df.ts <- na.omit(df.ts)
    df.ts$value <- .process.Date(df.ts$value)
    df.ts$yITL <- as.integer(factor(df.ts$individualID, levels=unique(df.ts$individualID)))
    names(df.ts) <- c("individualID", "type", "date", "yITL")


    ## GET THE DATA.FRAME FOR THE INDIVIDUAL LINES ##
    df.ind <- x@individuals[ordering,,drop=FALSE]

    ## keep only relevant individuals
    df.ind <- df.ind[unique(df.ts$individualID),,drop=FALSE]
    df.ind$individualID <- rownames(df.ind)
    df.ind$yITL <- 1:nrow(df.ind)


    ##df.full <- merge(df.ts, df.ind, by="individualID", all.x=TRUE, all.y=FALSE, sort=FALSE) # not needed


    ## BUILD THE PLOT ##
    ## base = time line
    out <- ggplot(df.ind)

    ## add horizontal lines for individuals
    if(is.null(colorBy)){
        out <- out + geom_hline(aes(yintercept=yITL),alpha=.3)
    } else{
        out <- out + geom_hline(aes_string(yintercept="yITL", colour=colorBy), alpha=I(.8), data=df.ind, show_guide=TRUE)
    }

    ## add points
    out <- out + geom_point(aes(x=date, y=yITL, shape=type), data=df.ts, ...)


    ## add indiv labels
    if(plotNames){
        out <- out + scale_y_discrete(name="Individuals", labels=df.ind$individualID)
    } else {
        out <- out + scale_y_discrete(name="Individuals",breaks=NULL)
    }

    ## THIS MAY NEED TESTING
    if(!is.null(periods)){
        if(class(periods)!="matrix"){
            warning("periods should be nx2 matrix of colnames")
        }
        else{
            for(i in 1:dim(periods)[1]){
                if(is.null(colorBy))
                    out <- out+geom_segment(aes_string(y="individualID",yend="individualID",x=periods[i,1],xend=periods[i,2]),size=1,lineend="round",alpha=.5)
                else
                    out <- out+geom_segment(aes_string(y="individualID",yend="individualID",x=periods[i,1],xend=periods[i,2],colour=colorBy),size=I(1),alpha=.5,lineend="round")

            }
        }
    }

    return(suppressWarnings(out))
} # end plotIndividualTimeline


