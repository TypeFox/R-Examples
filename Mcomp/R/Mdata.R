plot.Mdata <- function(x, xlim=c(tsp(x$x)[1],tsp(x$xx)[2]), ylim=range(x$x,x$xx),
                        main=x$sn, xlab="", ylab="", ...)
{
    freq = frequency(x$x)
    plot(ts(c(x$x, rep(NA, x$h)), end = tsp(x$x)[2] + x$h/freq, frequency = freq),ylim=ylim,
                xlim=xlim, ylab="", xlab="",...)
    if(nchar(x$description) > 70)
        title(main=list(main,cex=0.75,font=2))
    else
        title(main=list(main,cex=1,font=2))
    lines(ts(x$xx, start = tsp(x$x)[2] + 1/freq, frequency = freq),lt=1,col=2)
}

print.Mcomp <- function(x,...)
{
    n <- length(x)
    Type <- Period <- character(n)
    for(i in 1:n)
    {
        Type[i] <- x[[i]]$type
        if(Type[i]=="INDUSTRIAL")
            Type[i] <- "INDUST"
        if(strsplit(Type[i],"-")=="DEMOGRAPHI")
            Type[i] <- "DEMOGRAPHIC"
        Period[i] <- x[[i]]$period
    }
    if(length(unique(Type))==1)
        ty <- Type[1]
    if(length(unique(Period))==1)
        ty <- Period[1]

    if(length(unique(Type))>1 && length(unique(Period))>1)
    {
        cat(paste("M-Competition data:",n,"time series"),"\n\n")
        tab1 <- table(Period,Type)
        csum <- margin.table(tab1,2)
        tab2 <- as.table(rbind(tab1,Total=csum))
        rsum <- margin.table(tab2,1)
        tab3 <- as.table(cbind(tab2,Total=rsum),dnn=c("Period","Type of data"))
        names(dimnames(tab3))  <-  c("Period","Type of data")
        print(tab3)
    }
    else if(length(unique(Type))==1 && length(unique(Period))==1)
        print(table(Period,Type,dnn=c("Period","Type of data")))
    else
    {
        cat(paste("M-Competition data:",n,ty,"time series"),sep="","\n\n")
        print(table(Period,Type,dnn=c("Period","Type of data")))
    }
}

print.Mdata <- function(x,...)
{
    cat(paste(x$st),fill=TRUE,labels=paste("Series:"),sep="\n")
    cat(paste(x$type),fill=TRUE,labels=paste("Type of series:"),sep="\n")
    cat(paste(x$period),fill=TRUE,labels=paste("Period of series:"),sep="\n")
    cat(paste(x$description),fill=TRUE,labels=paste("Series description:"),sep="\n")
    cat("\nHISTORICAL data\n")
    print(x$x)
    cat("\nFUTURE data\n")
    print(x$xx)
}

Mcomp.sub <- function(x,getdata)
{
    n <- length(x)
    Type <- Period <- character(n)
    for(i in 1:n)
    {
        Type[i] <- x[[i]]$type
        if(Type[i]=="INDUSTRIAL")
            Type[i] <- "INDUST"
        if(strsplit(Type[i],"-")=="DEMOGRAPHI")
            Type[i] <- "DEMOGRAPHIC"
        Period[i] <- x[[i]]$period
    }
    getdatatable <- c("yearly","quarterly","monthly","other","macro","micro","industry",
                "finance","demographic","allother","macro1","macro2","micro1","micro2","micro3")
    if(is.character(getdata))
    {
        getdata <- getdatatable[charmatch(getdata,getdatatable)]
        if(length(getdata) != 1)
            stop("Ambiguous data type")
        else if(is.na(getdata))
            stop("Unknown data type")
    }

    if(getdata==1 | getdata=="yearly")
        choose <- (Period=="YEARLY")
    else if(getdata==4 | getdata=="quarterly")
        choose <- (Period=="QUARTERLY")
    else if(getdata==12 | getdata=="monthly")
        choose <- (Period=="MONTHLY")
    else if(getdata==111)
    {        
        j <- match(x111,names(x))
        choose <- rep(FALSE,length(x))
        choose[j] <- TRUE
    }
    else if(getdata=="other")
        choose <- 
        (Period=="OTHER")
    else if(getdata=="macro")
        choose <- (Type=="MACRO")
    else if(getdata=="micro")
        choose <- (Type=="MICRO")
    else if(getdata=="industry")
        choose <- (Type=="INDUSTRY" | Type=="INDUST")
    else if(getdata=="finance")
        choose <- (Type=="FINANCE")
    else if(getdata=="demographic")
        choose <- (Type=="DEMOGRAPHIC" | Type=="DEMOGR")
    else if(getdata=="allother")
        choose <- (Type=="OTHER")
    else if(getdata=="macro1")
        choose <- (Type=="MACRO1")
    else if(getdata=="macro2")
        choose <- (Type=="MACRO2")
    else if(getdata=="micro1")
        choose <- (Type=="MICRO1")
    else if(getdata=="micro2")
        choose <- (Type=="MICRO2")
    else if(getdata=="micro3")
        choose <- (Type=="MICRO3")
    else
        stop("Unknown type or period")
    if(sum(choose) == 0)
        stop("No data")

    return(x[choose])
}


subset.Mcomp <- function(x,cond1,cond2,...)
{
    if(is.character(cond1))
        cond1 <- tolower(cond1)
    M11=structure(Mcomp.sub(x,cond1),class="Mcomp")
    if(!missing(cond2))
    {
        if(is.character(cond2))
            cond2 <- tolower(cond2)
        M12=structure(Mcomp.sub(M11,cond2),class="Mcomp")
        return(M12)
    }
    else
        return(M11)
}
