#  *****************************************************************************
#   File : fdata.R
#         ************************************************************
#   Description : 
#       Functional Data class and methods.
#   Version : 2.0.1
#   Date : 2010-03-08
#         ************************************************************
#   Author : Julien Damon <julien.damon@gmail.com>
#   License : LGPL
#   URL: https://github.com/Looping027/far
#  *****************************************************************************

#  *****************************************************************************
#   Title : as.fdata
#         ************************************************************
#   Description : 
#       Generic method of 'fdata' constructor.
#       Select the right method using the class except for 
#       matrix and array.
#   Version : 2.0
#   Date : 2001-07-05
#  *****************************************************************************
as.fdata <- function(object,...)
{
    if (is.array(object))
    {
    if (is.matrix(object))
        return(as.fdata.matrix(object, ...))
    else
        return(as.fdata.array(object, ...))
    }
    else if (is.list(object) & (!is.data.frame(object)))
            return(as.fdata.list(object, ...))
        return(UseMethod("as.fdata"))
}

#  *****************************************************************************
#   Title : as.fdata.matrix
#         ************************************************************
#   Description : 
#       Method of 'fdata' constructor for 'matrix' object.
#   Version : 2.0
#   Date : 2003-06-02
#  *****************************************************************************
as.fdata.matrix <- function(object,...,col=0,p=0,dates,name=NULL)
{
    # if p=0  and col=0 the resulting object is almost the initial one 
    # with the class attribute
    if (p==0) 
        p <- nrow(object)
    if (col[1]==0)
    {
        if (is.null(colnames(object)))
        {
            if (missing(dates)) 
                colnames(object) <- 1:ncol(object)
            else 
                colnames(object) <- dates
        } else {
            if (!(missing(dates)))
                colnames(object) <- dates
        }
        if (is.null(rownames(object)))
            rownames(object) <- (0:(p-1))/p
        object <- list(object)
        if (is.null(name))
            names(object) <- "var"
        else
            names(object) <- name
        class(object) <- "fdata"
        return(object)
    } else {
        res <- list()
        if (is.null(name))
        {
            if ((is.null(colnames(object))) | (!(is.numeric(col))))
                colname <- col
            else
                colname <- colnames(object)[col]
        }
        else 
            colname <- name
        if (missing(dates))
            dates <- 1:(nrow(object)/p)
        if (is.null(rownames(object)))
            rowname <- (0:(p-1))/p
        else
            rowname <- rownames(object)[1:p]
        for (i in 1:length(col))
        {
            res[[i]] <- matrix(object[,col[i]],ncol=length(dates),nrow=p)
            dimnames(res[[i]]) <- list(rowname,dates)
        }
        names(res) <- colname
        class(res) <- "fdata"
        return(res)
    }
}

#  *****************************************************************************
#   Title : as.fdata.list
#         ************************************************************
#   Description : 
#       method of 'fdata' constructor for list
#       (just add the class attribute)
#   Version : 1.0
#   Date : 2003-06-02
#  *****************************************************************************
as.fdata.list <- function(object,...,dates,name=NULL)
{
    n <- length(object)
    for (i in 1:n)
    {
        p <- nrow(object[[i]])
        if (is.null(colnames(object[[i]])))
        {
            if (missing(dates)) 
                colnames(object[[i]]) <- 1:ncol(object[[i]])
            else 
                colnames(object[[i]]) <- dates
        } else {
            if (!(missing(dates)))
                colnames(object[[i]]) <- dates
        }
        if (is.null(rownames(object[[i]])))
            rownames(object[[i]]) <- (0:(p-1))/p
    }
    if (is.null(name)){
      if (is.null(names(object)))
        names(object) <- paste("var",1:n,sep="")
    }
    else
        names(object) <- name
    class(object) <- "fdata"
    return(object)
}

#  *****************************************************************************
#   Title : as.fdata.array
#         ************************************************************
#   Description : 
#       Method of 'fdata' constructor for 'array' object.
#   Version : 0.1
#   Date : 2010-03-08
#  *****************************************************************************
as.fdata.array <- function(object,...)
{
    # TODO : write a *real* version of that function
    as.fdata.matrix(as.matrix(object),...)
}
#  *****************************************************************************
#   Title : as.fdata.default
#         ************************************************************
#   Description : 
#       defualt method of 'fdata' constructor
#   Version : 2.0
#   Date : 2001-07-05
#  *****************************************************************************
as.fdata.default <- function(object,...)
{
    as.fdata.matrix(as.matrix(object),...)
}

#  *****************************************************************************
#   Title : print.fdata
#         ************************************************************
#   Description : 
#       print method for 'fdata' object
#   Version : 2.0
#   Date : 2001-07-05
#  *****************************************************************************
print.fdata <- function(x, ..., file="", append=TRUE)
{
    desc <- ncol(x[[1]])
    desc2 <- names(x)
    desc3 <- as.numeric(lapply(x,nrow))
    cat("Functionnal data\n\n",file=file,
        append=append)
    cat("Variable(s):",desc2,"\n",file=file,append=append)
    cat("Number of points per observation:",desc3,"\n",
        file=file,append=append)
    cat("Number of observations:",desc,"\n\n",file=file,
        append=append)
    invisible()
}

#  *****************************************************************************
#   Title : plot.fdata
#         ************************************************************
#   Description : 
#       plot method for 'fdata' object
#   Version : 1.1
#   Date : 2003-04-15
#  *****************************************************************************
plot.fdata <- function(x,...,date=1,xval=NULL,name=NULL,
                       main=NULL,whole=FALSE,separator=FALSE)
{
    n <- length(x)
    for (i in 1:n)
    {
        if (whole)
        {
            n2 <- ncol(x[[i]])
            n3 <- nrow(x[[i]])
            if (is.null(xval))
                x1 <- rep(as.numeric(rownames(x[[i]])),n2) +
                    rep(0:(n2-1),rep(n3,n2))
            else
                x1 <- xval
            y <- as.numeric(x[[i]])
            if (is.null(main)) {
              main1 <- paste(c(name,names(x)[i]),collapse="")
            } else {
              main1 <- main
            }
            plot(x=x1,y=y,type='l',xlab="time",
                ylab=names(x)[i],
                main=main1,...)
            if (separator)
            {
                min1 <- min(y)-10
                max1 <- max(y)+10
                for (j in 0:(n2-1))
                  { lines(x=rep(x1[1+j*n3],2),y=c(min1,max1),lty='dotted') }
            }            
        } else {
                if (is.null(xval))
                    x1 <- as.numeric(rownames(x[[i]]))
                else
                    x1 <- xval
        ndate <- length(date)
        for (j in 1:ndate)
        {
            if (is.null(main)) {
              main1 <- paste(c(name,names(x)[i]),collapse="")
            } else {
              main1 <- main
            }
            plot(x=x1,
                y=x[[i]][,date[j]],type='l',xlab="time",
                ylab=names(x)[i],
                main=main1,...)
        }
        }
    }
    invisible()   
}

#  *****************************************************************************
#   Title : summary.fdata
#         ************************************************************
#   Description : 
#       summary method for 'fdata' object
#   Version : 2.0
#   Date : 2001-07-06
#  *****************************************************************************
summary.fdata <- function(object, ...)
{
    # Generic function call in the differents cases
    summary.unifdata <- function(object)
    {
        res <- list()
        L1norm <- apply(abs(object),2,mean)
        L2norm <- sqrt(apply(object^2,2,mean))
        Linfnorm <- apply(abs(object),2,max)
        res$norm <- rbind(L1norm,L2norm,Linfnorm)
        dimnames(res$norm) <- list(c("L1 norm","L2 norm","Linf norm"),
                                    dimnames(object)[[2]])
        res$meannorm <- apply(res$norm,1,function(x)mean(x,na.rm=TRUE))
        return(res)
    }
    output <- list(NULL)
    length(output) <- length(object)
    names(output) <- names(object)
    for (i in 1:length(object))
    {
        output[[i]] <- summary.unifdata(object[[i]])
    }
    output$nbvar <- length(object)
    class(output) <- "summary.fdata"
    return(output)
}

#  *****************************************************************************
#   Title : print.summary.fdata
#         ************************************************************
#   Description : 
#       print method for 'summary.fdata' object
#   Version : 1.0
#   Date : 2001-03-16
#  *****************************************************************************
print.summary.fdata <- function(x, ..., file)
{
    if (missing(file))
        file <- NULL
    # Generic function call in the differents cases
    printuniv <- function(x,file=NULL)
    {
        if (is.null(file)) {
            cat("Mean of the norms:\n")
            print(x$meannorm)
        } else {
            cat("Mean of the norms:\n",file=file,append=TRUE)
            cat(names(x$meannorm),"\n",file=file,append=TRUE)
            cat(x$meannorm,"\n",file=file,append=TRUE)
        }
    return()
    }
    for (i in 1:as.numeric(x$nbvar))
    {
        if (is.null(file))
            cat("Variable: ",names(x)[i],"\n")
        else
            cat("Variable: ",names(x)[i],"\n",file=file,append=TRUE)
        printuniv(x[[i]],file)
        if (is.null(file))
            cat("\n")
        else
            cat("\n",file=file,append=TRUE)
    }
    invisible()
}

#  *****************************************************************************
#   Title : multplot.fdata
#         ************************************************************
#   Description : 
#       matplot method for 'fdata' object
#   Version : 1.0
#   Date : 2001-07-05
#  *****************************************************************************
multplot <- function (object, ...) UseMethod("multplot")

multplot.fdata <- function(object,date=1,xval=NULL,name=NULL,
        legend=FALSE,yleg,xlab=NULL,ylab=NULL,main=NULL,whole=FALSE,...)
{
    n <- length(object)
    n2 <- nrow(object[[1]])
    n3 <- ncol(object[[1]])
    yleglim <- missing(yleg)
    if (is.null(xlab)) xlab <- "time"
    if (is.null(ylab)) ylab <- ""
    if (whole) 
    {
        temp <- matrix(0,ncol=n,nrow=n3*n2)
        for (i in 1:n)
            temp[,i] <- as.numeric(object[[i]])
        n2 <- ncol(object[[1]])
        n3 <- nrow(object[[1]])
        if (is.null(xval))
            x <- rep(as.numeric(rownames(object[[1]])),n2) +
                rep(0:(n2-1),rep(n3,n2))
        else
            x <- xval
            if (is.null(main)) main2<-"" else main2<-main
        matplot(x=x,
            y=temp,type='l',...,
            xlab=xlab,
            ylab=ylab,
            main=main2,col=1:4,lty=c(1,2,4,5,6))
        if (legend)
            legend(x=min(x),y=max(temp),legend=names(object),col=1:4,
                   lty=c(1,2,4,5,6),bg="white")
    } else {
        temp <- matrix(0,ncol=n,nrow=n2)
        ndate <- length(date)
        if (is.null(xval))
            x <- as.numeric(rownames(object[[1]]))
        else
            x <- xval
        for (j in 1:ndate)
        {
            for (i in 1:n)
                temp[,i] <- (object[[i]])[,date[j]]
            if (is.null(main)) main2<-date[j] else main2<-main
            matplot(x=x,
                y=temp,type='l',xlab=xlab,
                ylab=ylab,
                main=main2,col=1:4,lty=c(1,2,4,5,6),...)
            if (yleglim) ylege <- max(temp)
            else ylege <- yleg
            if (legend)
                legend(x=min(x),y=ylege,legend=names(object),col=1:4,
                       lty=(c(1,2,4,5,6))[1:length(names(object))],bg="white")
        }
    }
    invisible()   
}


multplot.default <- function (object,...) 
{
    matplot(object,...)
}

#  *****************************************************************************
#   Title : fapply
#         ************************************************************
#   Description : 
#       Apply a function to each element of a fdata object
#   Version : 0.8
#   Date : 2002-12-17
#  *****************************************************************************
fapply <- function (data,FUN,row.names,...) 
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata"))
        stop("data is not of class fdata")
    FUN <- match.fun(FUN)
    f1 <- function(x) FUN(x, ...)
    if (missing(row.names)){
      f2 <- function(x) {
          nbcol <- ncol(x)
          res <- apply(x, 2, f1)
          dim(res) <- c(length(res)/nbcol, nbcol)        
          rownames(res) <- paste(1:nrow(res))
          colnames(res) <- colnames(x)
          return(res)
      }
    } else {
      f2 <- function(x) {
          nbcol <- ncol(x)
          res <- apply(x, 2, f1)
          dim(res) <- c(length(res)/nbcol, nbcol)
          rownames(res) <- row.names
          colnames(res) <- colnames(x)
          return(res)
      }
    }
    res <- lapply(data, f2)
    return(as.fdata(res, name = names(data)))
}

#  *****************************************************************************
#   Title : maxfdata
#         ************************************************************
#   Description : 
#       Compute the maximum of the day for each variable
#   Version : 2.0
#   Date : 2002-12-17
#  *****************************************************************************
maxfdata <- function(data)
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata"))
        stop("data is not of class fdata")
    fun1 <- function(x) max(x,na.rm=TRUE)
    return(fapply(data,fun1,row.names="max"))
}

#  *****************************************************************************
#   Title : is.na.fdata
#         ************************************************************
#   Description : 
#       Give the list of valid date
#   Version : 1.1
#   Date : 2007-04-07
#  *****************************************************************************
is.na.fdata <- function(x)
{
    n <- length(x)
    res <- matrix(1,ncol=ncol(x[[1]]),nrow=n)
    for (i in 1:n)
    {
        res[i,] <- c(apply(!is.na(x[[i]]),2,all))
    }
    colnames(res)<-colnames(x[[1]])
    rownames(res)<-names(x)
    return(res==0)
}

#  *****************************************************************************
#   Title : select.fdata
#         ************************************************************
#   Description : 
#       selection of fdata by dates
#   Version : 1.1
#   Date : 2003-05-15
#  *****************************************************************************
select.fdata <- function(data,date=NULL,name=NULL)
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata")) # test the type of data
        stop("data is not of class fdata")
    if (!is.null(name)) data <- data[name]
    if (!is.null(date)) data <- lapply(data,function(x,y=date)x[,y,drop=FALSE])
    class(data) <- "fdata"
    return(data)
}

#  *****************************************************************************
#   Title : date.fdata
#         ************************************************************
#   Description : 
#       gives the dates of a fdata
#   Version : 1.0
#   Date : 2001-07-09
#  *****************************************************************************
date.fdata <- function(data)
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata")) # test the type of data
        stop("data is not of class fdata")
    return(colnames(data[[1]]))
}

#  *****************************************************************************
#   Title : pred.persist
#         ************************************************************
#   Description : 
#       presistance forecasting for fdata
#   Version : 1.0
#   Date : 2001-07-09
#  *****************************************************************************
pred.persist <- function(data,x,na.rm=TRUE,
                label,positive=FALSE)
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata")) # test the type of data
        stop("data is not of class fdata")
    data <- data[x]
    class(data) <- "fdata"
    n <- ncol(data[[1]])
    p <- length(data)
    if (missing(label))
        label <- c(colnames(data[[1]])[-1],paste(n+1))
    else
        label <- c(colnames(data[[1]])[-1],label)
    for (i in 1:p)
        colnames(data[[i]]) <- label
    if (na.rm)
    {
        listobs <- c(apply(!is.na(data),2,all))        
        for (i in 1:p)
            data[[i]] <- (data[[i]])[,listobs]
    }
    if (positive)
        for (i in 1:p)
            data[[i]] <- (data[[i]]+abs(data[[i]]))/2
    return(data)
}
