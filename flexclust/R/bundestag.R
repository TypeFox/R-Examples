##
## Copyright (C) 2008 Friedrich Leisch
## $Id: bundestag.R 3 2013-06-12 10:06:43Z leisch $
##
bundestag <- function(year, second=TRUE, percent=TRUE, nazero=TRUE,
                      state=FALSE)
{
    year <- match.arg(as.character(year), c("2002", "2005", "2009"))

    getData <- function(year){
        tempenv <- new.env()
        load(system.file(paste("data/btw", year, ".RData", sep=""),
                         package="flexclust"),
             envir=tempenv)
        get(paste("btw", year, sep=""), envir=tempenv)
    }
        
    x <- getData(year)
    
    if(is.logical(state)){
        if(state) return(x$state)
    }
    else{
        y <- rep("other", nrow(x))
        ok <- grep(state, x$state)
        y[ok] <- as.character(x$state)[ok]
        return(as.factor(y))
    }

    if(second)
        p <- "2$"
    else
        p <- "1$"

    y <- x[,grep(p, colnames(x))]
    colnames(y) <- gsub(p, "", colnames(y))
    y <- as.matrix(y)
                        
    if(percent)
        y <- y/y[,"valid"]

    y <- y[,-grep("valid", colnames(y))]

    if(nazero) y[is.na(y)] <- 0

    y
}
