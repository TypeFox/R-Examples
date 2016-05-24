#############################################################################################
## File: PKfigure.R
## Author: Xiaoyong Sun
## Date: 10/12/2009
## Goal: PKreport
## Function list:
##    - PKreport.1: Exploratory data analysis
##    - PKreport.2: Individual plots
##    - PKreport.3: Goodness-of-fit plots
##    - PKreport.4: Structural model diagnostics
##    - PKreport.5: Residual model diagnostics
##    - PKreport.6: Parameters diagnostics
##    - PKreport.7: Covariate model
##    - PKreport.8: Random effects
##
#############################################################################################


################################################################################
## convert function: scatter plot, general plot -> xyplot or ggplot
##  - e.g. convert("xyplot", "qplot", lattice.list, ggplot.list, plot.dir, pdata, paste.string-for complex string)
################################################################################

## INPUT
##    - lattice.call: xyplot
##    - lattice.string: command string without lattice.call(from paste)
##    - ggplot.call: ggplot
##    - ggplot.string: command string without lattice.call(from paste)
##    - plot.dir: folder to save plot
##    - paste.string: for those not generated from addList(), complex command,
##                    command string with lattice/ggplot call
## OUTPUT: figures
generate.plot <- function(lattice.call, lattice.string=NULL, ggplot.call, 
                  ggplot.string=NULL, plot.dir, pdata, paste.string=TRUE)
{
    general.list = .pkplot$getGlobalConfig()
    save.call <- general.list$save.format

    ## save code scripts
    if (.pkplot$getConfig("package")==0 || .pkplot$getConfig("package")==2)
    {
        if (!is.null(lattice.string))
        {
            if (paste.string)
            {
                lattice.code <- paste(lattice.call, "(", lattice.string, ")", sep="")
            }
            else
            {
                lattice.code <- lattice.string
            }
            tmp.lattice <- list()
            tmp.lattice[[1]] <- lattice.code
            tmp.lattice[[2]] <- plot.dir
            .pkplot$setPKCode(tmp.lattice)

            try1 <- try(eval(parse(text=lattice.code)))
            if (!inherits(try1, "try-error"))
            {
                for (i in 1:length(save.call))
                {
                    fig.name <- addFigName(plot.dir, save.call[i])
                    do.call(save.call[i], list(file=paste("lattice", fig.name, sep="_"), width=general.list$width,
                              height=general.list$height))
                    print(eval(parse(text=lattice.code)))
                    dev.off()
                }
            }
        }
    }

    ## save code scripts
    if (.pkplot$getConfig("package")==1 || .pkplot$getConfig("package")==2)
    {
        if (!is.null(ggplot.string))
        {
            if (paste.string)
            {
                ggplot.code <- paste(ggplot.call, "(", ggplot.string, ")", sep="")
            }
            else
            {
                ggplot.code <- ggplot.string
            }

            tmp.lattice <- list()
            tmp.lattice[[1]] <- ggplot.code
            tmp.lattice[[2]] <- plot.dir
            .pkplot$setPKCode(tmp.lattice)

            try1 <- try(eval(parse(text=ggplot.code)))
            if (!inherits(try1, "try-error"))
            {
                for (i in 1:length(save.call))
                { 
                  fig.name <- addFigName(plot.dir, save.call[i])           
                  do.call(save.call[i], list(file=paste("ggplot", fig.name, sep="_"), width=general.list$width,
                            height=general.list$height))
                  print(eval(parse(text=ggplot.code)))
                  dev.off()
                }
            }
        }
    }

    invisible(NULL)
}


addList1 <- function(lattice.list, ggplot.list, x)
{
        lattice.list$x <- formula(paste("~", x, sep=""))

        ggplot.list$x <- as.name(x)
        
        lattice.string <- paste("xlab=\"", x, "\"", sep="")
        sapply(1:length(lattice.list), function(i)
              {
                  if (!is.character(lattice.list[[i]]))
                  {
                      lattice.string <<- paste(lattice.string, ",",
                                         names(lattice.list[i]), "=", lattice.list[i], sep="")
                  }
                  else
                  {
                      lattice.string <<- paste(lattice.string, ",",
                                         names(lattice.list[i]), "=\"", lattice.list[i], "\"",sep="")
                  }
              })

        ggplot.string <- paste("xlab=\"", x, "\"", sep="")
        sapply(1:length(ggplot.list), function(i)
              {
                  if (!is.character(ggplot.list[[i]]))
                  {
                      ggplot.string <<- paste(ggplot.string, ",",
                                         names(ggplot.list[i]), "=", ggplot.list[i], sep="")
                  }
                  else
                  {
                      ggplot.string <<- paste(ggplot.string, ",",
                                         names(ggplot.list[i]), "=\"", ggplot.list[i], "\"",sep="")
                  }
              })

        return(list(lattice.string=lattice.string, ggplot.string=ggplot.string))
}
addList1.old <- function(lattice.list, ggplot.list, x)
{
        lattice.list$x <- formula(paste("~", x, sep=""))
        lattice.list$xlab <- x

        ggplot.list$x <- as.name(x)
        ggplot.list$xlab <- x

        return(list(lattice.list=lattice.list, ggplot.list=ggplot.list))
}

addList2 <- function(lattice.list, ggplot.list, x, y, cond)
{
        if (!is.null(x) && !is.null(y))
        {
            if (missing(cond))  lattice.list$x <- formula(paste(y, "~", x, sep=""))
            else lattice.list$x <- formula(paste(y, "~", x, "|", cond, sep=""))
        }
        
        tmp.string <- paste(names(lattice.list), lattice.list, sep="=")
        lattice.string <- paste(tmp.string, collapse=",")
        
        if (!is.null(x) && !is.null(y))
        {
            lattice.string <- paste("xlab=\"", x, "\"", sep="")
            lattice.string <- paste(lattice.string, ",ylab=\"", y, "\"", sep="")
        }
        else
        {
            lattice.string <- paste("xlab=\"", NULL, "\"", sep="")
            lattice.string <- paste(lattice.string, ",ylab=\"", NULL, "\"", sep="")
        }
        sapply(1:length(lattice.list), function(i)
              {
                  if (!is.character(lattice.list[[i]]) || length(lattice.list[[i]]) >=2 )
                  {
                      lattice.string <<- paste(lattice.string, ",",
                                         names(lattice.list[i]), "=", lattice.list[i], sep="")
                  }
                  else
                  {
                      lattice.string <<- paste(lattice.string, ",",
                                         names(lattice.list[i]), "=\"", lattice.list[i], "\"",sep="")
                  }
              })

        if (!is.null(x) && !is.null(y))
        {
            ggplot.list$x <- as.name(x)
            ggplot.list$y <- as.name(y)
            ggplot.string <- paste("xlab=\"", x, "\"", sep="")
            ggplot.string <- paste(ggplot.string, ",ylab=\"", y, "\"", sep="")
        }
        else
        {
            ## start ggplot.list
            ggplot.string <- paste("xlab=\"", NULL, "\"", sep="")
            ggplot.string <- paste(ggplot.string, ",ylab=\"", NULL, "\"", sep="")
        }

        ggplot.list$se <- FALSE

        sapply(1:length(ggplot.list), function(i)
              {
                  if (!is.character(ggplot.list[[i]])|| length(lattice.list[[i]]) >=2 )
                  {
                      ggplot.string <<- paste(ggplot.string, ",",
                                         names(ggplot.list[i]), "=", ggplot.list[i], sep="")
                  }
                  else
                  {
                      ggplot.string <<- paste(ggplot.string, ",",
                                         names(ggplot.list[i]), "=\"", ggplot.list[i], "\"",sep="")
                  }
              })

        return(list(lattice.string=lattice.string, ggplot.string=ggplot.string))
}

## Create fig name for lattice and ggplot sequentially
#   input: - plot.dir
#          - save.format
#          - i: fig no
#   output: a fig name for figure comprising the fig no from cache
addFigName <- function(plot.dir, save.format, i=0)
{
        if (i==0) i <- .pkplot$getPKCodeLen()

        fig.name <- paste(i, plot.dir, sep="_")
        if(save.format=="win.metafile")
        {
          fig.name <- paste(fig.name, "emf", sep=".")
        }
        else
        {
          fig.name <- paste(fig.name, save.format, sep=".")
        }
        return(fig.name)
}

addCodeNote <- function(newstr)
{
        ## add word explain for R code
     word.repeat <- .pkplot$getConfig("package")
     if (word.repeat == 2)
     {
        .pkplot$setPKCodeNote(paste(newstr, " (lattice)", sep=""))
        .pkplot$setPKCodeNote(paste(newstr, " (ggplot2)", sep=""))
     }
     else
     {
        .pkplot$setPKCodeNote(newstr)
     }

}

## Goal: to make term1 and term2 have equal length
## example: PRED <- PRED, RES <- c(CRES, RES)
## Note: additional part is from first [1] repeat
twoTermEqualLen <- function(myterm1, myterm2)
{
    if (length(myterm1)==length(myterm2)) return(list(term1=myterm1, term2=myterm2))
    if (length(myterm1) > length(myterm2))
    {
        len.diff <- length(myterm1) - length(myterm2)
        myterm2 <- c(myterm2, rep(myterm2[1], len.diff))
    }
    else
    {
        len.diff <- length(myterm2) - length(myterm1)
        myterm1 <- c(myterm1, rep(myterm1[1], len.diff))
    }
    return(list(term1=myterm1, term2=myterm2))
}

################################################################################
## PK.run.1: simple EDA
################################################################################
PKreport.1 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")

    ## global config
    match.term <- .pkplot$getTerm()
    general.list <- .pkplot$getGlobalConfig()
    lattice.list <- .pkplot$getHistGraph("lattice")
    ggplot.list <- .pkplot$getHistGraph("ggplot")

    #####################
    ## univariate graph - histogram
    #####################
    plot.dir <- general.list$univar.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup
    lattice.call <- "histogram"
    ggplot.call <- "qplot"
    
    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)

    for (i in 1:length(colnames(pdata)))
    {

        newstr <- paste("Distribution of ", colnames(pdata)[i], sep="")
        addCodeNote(newstr)
        ggplot.list.tmp <- ggplot.list

        bin.tmp <- round(diff(range(pdata[,colnames(pdata)[i]]))/30,5)
        if (bin.tmp != 0) ggplot.list.tmp$binwidth <- bin.tmp
        else ggplot.list.tmp$binwidth <- diff(range(pdata[,colnames(pdata)[i]]))/30

        my.list <- addList1(lattice.list, ggplot.list.tmp, colnames(pdata)[i])

        generate.plot(lattice.call, my.list$lattice.string, 
                      ggplot.call, my.list$ggplot.string, 
                      plot.dir, pdata)
    }
    
    setwd(tmp.dir)

    #####################
    ## bivariate graph
    #####################

    plot.dir <- general.list$bivar.dir
    if(!is.null(plot.dir)&& (!file.exists(plot.dir))) dir.create(plot.dir)
    tmp.dir <- getwd()
    setwd(paste(tmp.dir, plot.dir, sep="/"))
    
    lattice.call <- "xyplot"

    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")
    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)
    
    ## get rid of se band in ggplot2
    ggplot.list$se <- FALSE
    
    # only do bivariate figure for terms
    match.term <- .pkplot$getTerm()
    bipara <- unique(unlist(match.term))
    for (i in 1:(length(bipara)-1))  #x
    {
        for (j in 2:(length(bipara))) #y
        {
          newstr <- paste(bipara[j], " vs ", bipara[i], sep="")
          addCodeNote(newstr)
          my.list <- addList2(lattice.list, ggplot.list, bipara[j], bipara[i])

          generate.plot(lattice.call, my.list$lattice.string, 
                        ggplot.call, my.list$ggplot.string, 
                        plot.dir, pdata)
        }
    }
    invisible(NULL)
}
################################################################################
################################################################################
## Individual plots
PKreport.2 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")

    ## require terms:
    match.term <- .pkplot$getTerm()
    require.term <- c("DV", "PRED", "IPRE", "IDV")
    PK.match <- match(require.term, names(match.term))
    if(length(PK.match[is.na(PK.match)]) > 0)
        stop(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined in var.name!", sep=""))

    ## global config
    general.list <- .pkplot$getGlobalConfig()
    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")

    plot.dir <- general.list$ind.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup
    lattice.call <- "xyplot"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)

    # TODO
    # Individual plots (IDV vs DV/PRED/IPRE | ID)
    if (  (length(.pkplot$getTerm()$PRED)!=0) && (length(.pkplot$getTerm()$IDV)!=0) &&
            (length(.pkplot$getTerm()$IPRE)!=0) )
    {
        for (i in 1: length(match.term$IDV))
        {
            newstr <- paste("individual plots for ", match.term$IDV[i], " vs ", match.term$DV,
                            "/", match.term$PRED,
                            "/", match.term$IPRE, " | ", match.term$ID, sep="")
            addCodeNote(newstr)

            ggplot.layout <- .pkplot$getScatterGraph("others")$ind.layout
            
            tmp.id <- unique(pdata[,match.term$ID])            
            if (length(tmp.id)>25) 
            {
                tmp.index <- which(pdata[,match.term$ID] %in% tmp.id[1:25])
                tmp.pdata <- pdata[tmp.index,] 
            }
            else
            {
                tmp.pdata <- pdata
            }
 
            tmp.string1 <- paste("new.data <- pdata[c(\"", match.term$ID,
                            "\", \"", match.term$IDV,
                            "\", \"", match.term$DV,
                            "\", \"", match.term$PRED,
                            "\", \"", match.term$IPRE, "\")]", sep="")
            tmp.string2 <- paste("new.data2 <- data.frame(IDV=new.data[, \"", match.term$IDV, "\"], ID=new.data[,\"",
                            match.term$ID, "\"], y=c(new.data[,\"",
                            match.term$DV, "\"], new.data[,\"", match.term$IPRE,
                            "\"], new.data[,\"", match.term$IPRE,
                            "\"]), label=c(rep(\"", match.term$DV,
                            "\", nrow(new.data)), rep(\"", match.term$PRED,
                            "\", nrow(new.data)), rep(\"", match.term$IPRE,
                            "\", nrow(new.data))))", sep="")
            tmp.string3.lattice <- paste("xyplot(y~IDV|factor(ID), group=label, data=new.data2,
                                     auto.key=list(space=\"right\"), as.table=TRUE, layout=c(",
                                     ggplot.layout[1], ",", ggplot.layout[2],
                                     "), xlab=\"IDV\", ylab=\"DV/PRED/IPRE\")", sep="")
            tmp.string3.ggplot <- paste("qplot(x=IDV, y=y, colour=label, data=new.data2,xlab=\"IDV\",
                                  ylab=\"DV/PRED/IPRE\")+ facet_wrap(~ID, ncol =",
                                  as.numeric(ggplot.layout[1]), ")", sep="")

            lattice.string <- paste(tmp.string1, tmp.string2, tmp.string3.lattice, sep=";")
            ggplot.string <- paste(tmp.string1, tmp.string2, tmp.string3.ggplot, sep=";")
            generate.plot(lattice.call, lattice.string,
                        ggplot.call, ggplot.string, plot.dir, tmp.pdata, paste.string=FALSE)

        }
    }
}
################################################################################
################################################################################
## Goodness-of-fit plots
PKreport.3 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")

    ## require terms:
    match.term <- .pkplot$getTerm()
    require.term <- c("DV", "PRED", "IPRE", "IDV")
    PK.match <- match(require.term, names(match.term))
    if(length(PK.match[is.na(PK.match)]) > 0)
        stop(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined in var.name!", sep=""))

    ## global config
    general.list <- .pkplot$getGlobalConfig()
    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")

    plot.dir <- general.list$gof.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))
    
    ## global setup
    lattice.call <- "xyplot"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)

    # DV vs PRED
    if (  (length(.pkplot$getTerm()$DV)!=0) && (length(.pkplot$getTerm()$PRED)!=0) )
    {
          newstr <- paste(match.term$DV, "vs", match.term$PRED, sep=" ")
          addCodeNote(newstr)

          my.list <- addList2(lattice.list, ggplot.list, match.term$DV, match.term$PRED)
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, 
          plot.dir, pdata)

    }
    # DV vs IPRED
    if (  (length(.pkplot$getTerm()$DV)!=0) && (length(.pkplot$getTerm()$IPRE)!=0) )
    {
          newstr <- paste(match.term$DV, "vs", match.term$IPRE, sep=" ")
          addCodeNote(newstr)
          my.list <- addList2(lattice.list, ggplot.list, match.term$DV, match.term$IPRE)
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, 
          plot.dir, pdata)

    }

    # DV vs Predictions(PRED / IPRED, x-axis) (dv.vs.pred.ipred).
    if (  (length(.pkplot$getTerm()$DV)!=0) && (length(.pkplot$getTerm()$IPRE)!=0)
                && (length(.pkplot$getTerm()$PRED)!=0))
    {
            newstr <- paste(match.term$DV, " vs ", match.term$PRED, "/", match.term$IPRE, sep="")
            addCodeNote(newstr)

            tmp.string1 <- paste("new.data <- pdata[,c(\"", match.term$DV, "\",\"", match.term$PRED, "\",\"",
                                 match.term$IPRE, "\")]", sep="")
            tmp.string2 <- paste("new.data2 <- data.frame(DV=new.data[,\"", match.term$DV, "\"], y=c(new.data[,\"",
                                match.term$PRED, "\"], new.data[,\"", match.term$IPRE,
                                "\"]), label=c(rep(\"", match.term$PRED, "\", nrow(new.data)), rep(\"",
                                match.term$IPRE, "\", nrow(new.data))))", sep="")
            tmp.string3.lattice <- "xyplot(y~DV, group=label, data=new.data2, xlab=\"DV\", ylab=\"PRED/IPRE\")"
            tmp.string3.ggplot <- "qplot(x=DV, y=y, colour=label, data=new.data2, xlab=\"DV\", ylab=\"PRED/IPRE\")"

        my.list$lattice.string <- paste(tmp.string1, tmp.string2, tmp.string3.lattice, sep=";")
        my.list$ggplot.string <- paste(tmp.string1, tmp.string2, tmp.string3.ggplot, sep=";")
        generate.plot(lattice.call, my.list$lattice.string, 
        ggplot.call, my.list$ggplot.string, plot.dir, 
        pdata, paste.string=FALSE)
    }

    # PRED vs IDV
    if (  (length(.pkplot$getTerm()$PRED)!=0) && (length(.pkplot$getTerm()$IDV)!=0) )
    {

          for (i in 1: length(match.term$IDV))
          {
              newstr <- paste(match.term$PRED, "vs", match.term$IDV[i], sep=" ")
              addCodeNote(newstr)
              my.list <- addList2(lattice.list, ggplot.list, match.term$IDV[i], match.term$PRED)
              generate.plot(lattice.call, my.list$lattice.string, 
              ggplot.call, my.list$ggplot.string, plot.dir, pdata)

          }

    }
    
    # IPRE vs IDV
    if (  (length(.pkplot$getTerm()$IPRE)!=0) && (length(.pkplot$getTerm()$IDV)!=0) )
    {
          for (i in 1: length(match.term$IDV))
          {
              newstr <- paste(match.term$IPRE, "vs", match.term$IDV[i], sep=" ")
              addCodeNote(newstr)
              my.list <- addList2(lattice.list, ggplot.list, match.term$IDV[i], match.term$IPRE)
              generate.plot(lattice.call, my.list$lattice.string, 
              ggplot.call, my.list$ggplot.string, plot.dir, pdata)
          }

    }

    # Predictions vs IDV (PRED / IPRED, y-axis)
    if (  (length(.pkplot$getTerm()$IDV)!=0) && (length(.pkplot$getTerm()$IPRE)!=0)
                && (length(.pkplot$getTerm()$PRED)!=0))
    {

            newstr <- paste(match.term$IDV, " vs ", match.term$PRED, "/", match.term$IPRE, sep="")
            addCodeNote(newstr)

            tmp.string1 <- paste("new.data <- pdata[,c(\"", match.term$IDV, "\",\"", match.term$PRED, "\",\"",
                                 match.term$IPRE, "\")]", sep="")
            tmp.string2 <- paste("new.data2 <- data.frame(IDV=new.data[,\"", match.term$IDV, "\"], y=c(new.data[,\"",
                                match.term$PRED, "\"], new.data[,\"", match.term$IPRE,
                                "\"]), label=c(rep(\"", match.term$PRED, "\", nrow(new.data)), rep(\"",
                                match.term$IPRE, "\", nrow(new.data))))", sep="")
            tmp.string3.lattice <- "xyplot(y~IDV, group=label, data=new.data2, xlab=\"IDV\", ylab=\"PRED/IPRE\")"
            tmp.string3.ggplot <- "qplot(x=IDV, y=y, colour=label, data=new.data2, xlab=\"IDV\", ylab=\"PRED/IPRE\")"

            my.list$lattice.string <- paste(tmp.string1, tmp.string2, tmp.string3.lattice, sep=";")
            my.list$ggplot.string <- paste(tmp.string1, tmp.string2, tmp.string3.ggplot, sep=";")
            generate.plot(lattice.call, my.list$lattice.string, 
            ggplot.call, my.list$ggplot.string, plot.dir, pdata, paste.string=FALSE)

    }

}

# Structural model diagnostics
PKreport.4 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")

    ## require terms:
    match.term <- .pkplot$getTerm()
    require.term <- c("DV", "PRED", "IPRE", "IDV", "WRES", "COV")
    PK.match <- match(require.term, names(match.term))
    if(length(PK.match[is.na(PK.match)]) > 0)
        stop(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined in var.name!", sep=""))

    ## global config
    general.list <- .pkplot$getGlobalConfig()
    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")
    
    ## create folder now
    plot.dir <- general.list$struct.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup
    lattice.call <- "xyplot"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)
    ggplot.layout <- .pkplot$getScatterGraph("others")$ind.layout

    #PRED vs DV|IDV:
    if (  (length(.pkplot$getTerm()$PRED)!=0) && (length(.pkplot$getTerm()$IDV)!=0) )
    {
        for (i in 1: length(match.term$IDV))
        {
            newstr <- paste(match.term$PRED," vs ", match.term$DV, " | ", match.term$IDV[i], sep="")
            addCodeNote(newstr)
            
            tmp.id <- unique(pdata[,match.term$IDV[i]])            
            if (length(tmp.id)>25) 
            {
                tmp.index <- which(pdata[,match.term$IDV[i]] %in% tmp.id[1:25])
                tmp.pdata <- pdata[tmp.index,] 
            }
            else
            {
                tmp.pdata <- pdata
            }

            lattice.string <- paste("xyplot(", match.term$PRED, "~", match.term$DV, "|factor(", match.term$IDV[i],
                         "), layout=c(", ggplot.layout[1], ",", ggplot.layout[2], "), data=pdata, as.table=TRUE)",
                         sep="")
            tmp.string1 <- paste("pdata[,\"", match.term$IDV[i], "\"] <- factor(pdata[,\"", match.term$IDV[i], "\"])", sep="")
            tmp.string2 <- paste("ggplot(pdata)+geom_point(aes(", match.term$DV, ",", match.term$PRED, "))",
                          "+facet_wrap(~", match.term$IDV[i], ", ncol=", as.numeric(ggplot.layout[1]),")",
                          sep="")
            ggplot.string <- paste(tmp.string1, tmp.string2, sep=";")
            generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                          plot.dir, tmp.pdata, paste.string=FALSE)

        }
    }

    #IPRED vs DV|IDV
    if (  (length(.pkplot$getTerm()$IPRE)!=0) && (length(.pkplot$getTerm()$IDV)!=0) )
    {
        for (i in 1: length(match.term$IDV))
        {
            newstr <- paste("IPRE vs DV | ", match.term$IDV[i], sep="")
            addCodeNote(newstr)

            tmp.id <- unique(pdata[,match.term$IDV[i]])            
            if (length(tmp.id)>25) 
            {
                tmp.index <- which(pdata[,match.term$IDV[i]] %in% tmp.id[1:25])
                tmp.pdata <- pdata[tmp.index,] 
            }
            else
            {
                tmp.pdata <- pdata
            }

            lattice.string <- paste("xyplot(", match.term$IPRE, "~", match.term$DV, "|factor(", match.term$IDV[i],
                         "), layout=c(", ggplot.layout[1], ",", ggplot.layout[2], "), data=pdata, as.table=TRUE)",
                         sep="")
            tmp.string1 <- paste("pdata[,\"", match.term$IDV[i], "\"] <- factor(pdata[,\"", match.term$IDV[i], "\"])", sep="")
            tmp.string2 <- paste("ggplot(pdata)+geom_point(aes(", match.term$DV, ",", match.term$IPRE, "))",
                          "+facet_wrap(~", match.term$IDV[i], ", ncol=", as.numeric(ggplot.layout[1]),")",
                          sep="")
            ggplot.string <- paste(tmp.string1, tmp.string2, sep=";")
            generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                          plot.dir, tmp.pdata, paste.string=FALSE)
        }
    }

    #WRES vs IDV:
    if (  (length(.pkplot$getTerm()$WRES)!=0) && (length(.pkplot$getTerm()$IDV)!=0) )
    {
      for (i in 1: length(match.term$IDV))
      {
          newstr <- paste(match.term$WRES, " vs ", match.term$IDV[i], sep="")
          addCodeNote(newstr)
          my.list <- addList2(lattice.list, ggplot.list, match.term$IDV[i], match.term$WRES)
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, plot.dir, pdata)

      }
    }

    #WRES vs PRED:
    if (  (length(.pkplot$getTerm()$WRES)!=0) && (length(.pkplot$getTerm()$PRED)!=0) )
    {    
      newstr <- paste(match.term$WRES, " vs ", match.term$PRED, sep="")
      addCodeNote(newstr)
      my.list <- addList2(lattice.list, ggplot.list, match.term$PRED, match.term$WRES)
      generate.plot(lattice.call, my.list$lattice.string, 
      ggplot.call, my.list$ggplot.string, plot.dir, pdata)
    
    ##WRES vs PRED (bw):   
      #newstr <- paste(match.term$WRES, " vs ", match.term$PRED, "_bw", sep="")
      #addCodeNote(newstr)  
      #tmp.string1 <- paste("pdata[,\"", match.term$PRED, "\"] <- factor(pdata[,\"",
      #                    match.term$IDV, "\"])", sep="")
      #tmp.string2 <- paste("bwplot(", match.term$WRES, "~", match.term$PRED, ",data=",
      #                   quote(pdata), ",xlab=\"", match.term$PRED, "\",ylab=\"",
      #                   match.term$WRES, "\")", sep="")
      #lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
  
      #tmp1 <- paste("ggplot(pdata, aes(factor(", match.term$PRED, "),", match.term$WRES, sep="")
      #tmp2 <- paste(")) + geom_boxplot() +", "labs(x=\"", match.term$PRED, "\", y=\"", match.term$WRES, "\")", sep="")
      #tmp3 <- paste(tmp1, tmp2, sep="")
      #ggplot.string <- paste(tmp.string1, tmp3, sep=";")
      #
      #generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string, plot.dir, FALSE)
    
    }
    
    #PRED vs DV|Covariates
    if (  (length(.pkplot$getTerm()$DV)!=0) && (length(.pkplot$getTerm()$PRED)!=0) &&
          (length(.pkplot$getTerm()$COV)!=0)  )
    {     
      for (i in 1: length(match.term$COV))
      {
          newstr <- paste(match.term$PRED, " vs ", match.term$DV, " | ", match.term$COV[i], sep="")
          addCodeNote(newstr)
            
            tmp.id <- unique(pdata[,match.term$COV[i]])            
            if (length(tmp.id)>25) 
            {
                tmp.index <- which(pdata[,match.term$COV[i]] %in% tmp.id[1:25])
                tmp.pdata <- pdata[tmp.index,] 
            }
            else
            {
                tmp.pdata <- pdata
            }

          lattice.string <- paste("xyplot(", match.term$PRED, "~", match.term$DV, "|factor(", match.term$COV[i],
                         "), layout=c(", ggplot.layout[1], ",", ggplot.layout[2], "), data=pdata, as.table=TRUE)",
                         sep="")
          tmp.string1 <- paste("pdata[,\"", match.term$COV[i], "\"] <- factor(pdata[,\"", match.term$COV[i], "\"])", sep="")
          tmp.string2 <- paste("ggplot(pdata)+geom_point(aes(", match.term$DV, ",", match.term$PRED, "))",
                          "+facet_wrap(~", match.term$COV[i], ", ncol=", as.numeric(ggplot.layout[1]),")",
                          sep="")
          ggplot.string <- paste(tmp.string1, tmp.string2, sep=";")
          generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                          plot.dir, tmp.pdata, paste.string=FALSE)
  
      }
    }
    #IPRED vs DV|Covariates dependent variable (DV, y-axis) plotted against individual predictions (IPRED, x-axis), conditioned on covariates (dv.vs.ipred.by.cov).
    if (  (length(.pkplot$getTerm()$DV)!=0) && (length(.pkplot$getTerm()$IPRE)!=0) &&
          (length(.pkplot$getTerm()$COV)!=0)  )
    {      
      for (i in 1: length(match.term$COV))
      {
          newstr <- paste(match.term$IPRE, " vs ", match.term$DV, " | ", match.term$COV[i], sep="")
          addCodeNote(newstr)
          
            tmp.id <- unique(pdata[,match.term$COV[i]])            
            if (length(tmp.id)>25) 
            {
                tmp.index <- which(pdata[,match.term$COV[i]] %in% tmp.id[1:25])
                tmp.pdata <- pdata[tmp.index,] 
            }
            else
            {
                tmp.pdata <- pdata
            }

          lattice.string <- paste("xyplot(", match.term$IPRE, "~", match.term$DV, "|factor(", match.term$COV[i],
                         "), layout=c(", ggplot.layout[1], ",", ggplot.layout[2], "), data=pdata, as.table=TRUE)",
                         sep="")
          tmp.string1 <- paste("pdata[,\"", match.term$COV[i], "\"] <- factor(pdata[,\"", match.term$COV[i], "\"])", sep="")
          tmp.string2 <- paste("ggplot(pdata)+geom_point(aes(", match.term$DV, ",", match.term$IPRE, "))",
                          "+facet_wrap(~", match.term$COV[i], ", ncol=", as.numeric(ggplot.layout[1]),")",
                          sep="")
          ggplot.string <- paste(tmp.string1, tmp.string2, sep=";")
          generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                          plot.dir, tmp.pdata, paste.string=FALSE)
      }
    }

    ## global boxplot setup
    lattice.call <- "bwplot"
    ggplot.call <- "qplot"
    
    lattice.list <- list(data=quote(pdata))
    ggplot.list <- list(data=quote(pdata), geom=c("boxplot"))
    
    #WRES vs IDV (boxplot)
    if (  (length(.pkplot$getTerm()$WRES)!=0) && (length(.pkplot$getTerm()$DV)!=0) )
    {       
      for (i in 1: length(match.term$IDV))
      {
          newstr <- paste(match.term$WRES, " vs ", match.term$IDV[i], "_bw", sep="")
          addCodeNote(newstr)
          tmp.string1 <- paste("pdata[,\"", match.term$IDV[i], "\"] <- factor(pdata[,\"",
                          match.term$IDV, "\"])", sep="")
          tmp.string2 <- paste("bwplot(", match.term$WRES, "~ factor(", match.term$IDV[i], "),data=",
                         quote(pdata), ",xlab=\"", match.term$IDV[i], "\",ylab=\"",
                         match.term$WRES, "\")", sep="")
          lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
  
          tmp1 <- paste("ggplot(pdata, aes(factor(", match.term$IDV[i], "),", match.term$WRES, sep="")
          tmp2 <- paste(")) + geom_boxplot() +", "labs(x=\"", match.term$IDV[i], "\", y=\"", match.term$WRES, "\")", sep="")
          tmp3 <- paste(tmp1, tmp2, sep="")
          ggplot.string <- paste(tmp.string1, tmp3, sep=";")
  
          generate.plot(lattice.call, lattice.string, 
          ggplot.call, ggplot.string, plot.dir, tmp.pdata, FALSE)
  
      }
     }
}

# Residual model diagnostics
PKreport.5 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")
    ## require terms:
    match.term <- .pkplot$getTerm()
    require.term <- c("PRED", "WRES", "COV")
    PK.match <- match(require.term, names(match.term))
    if(length(PK.match[is.na(PK.match)]) > 0)
        stop(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined in var.name!", sep=""))

    .pkplot$setPKData(pdata)
    ## global config
    general.list <- .pkplot$getGlobalConfig()
    
    ##create folder
    plot.dir <- general.list$resid.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup for hist
    lattice.list <- .pkplot$getHistGraph("lattice")
    ggplot.list <- .pkplot$getHistGraph("ggplot")

    lattice.call <- "histogram"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)

    if (  length(.pkplot$getTerm()$WRES)!=0 && length(.pkplot$getTerm()$ID)!=0 )
    { 
      # Distribution of WRES (hist)
      newstr <- paste("Distribution of ", match.term$WRES, sep="")
      addCodeNote(newstr)
      ggplot.list.tmp <- ggplot.list
      
      bin.tmp <- round(diff(range(pdata[,match.term$WRES]))/30,5)
      if (bin.tmp != 0) ggplot.list.tmp$binwidth <- bin.tmp
      else  ggplot.list.tmp$binwidth <- diff(range(pdata[,match.term$WRES]))/30
      
      my.list <- addList1(lattice.list, ggplot.list.tmp, match.term$WRES)
      generate.plot(lattice.call, my.list$lattice.string, 
      ggplot.call, my.list$ggplot.string, plot.dir, pdata)

      # Distribution of WRES (QQ)
      newstr <- paste("Distribution of ", match.term$WRES, "(QQ)", sep="")
      addCodeNote(newstr)
      my.list$lattice.string <- paste("qqmath(~", match.term$WRES, ",data=", quote(pdata), 
             ")",
            #", prepanel = prepanel.qqmathline, panel = function(x, ...) {panel.qqmathline(x, ...);panel.qqmath(x, ...)})",
             sep="")
      #my.list$ggplot.string <- paste("ggplot(", quote(pdata), ") + geom_point(aes(sample =",
      #                              match.term$WRES, "), stat = \"qq\", quantiles = ppoints(100))",
      #                              sep="")
      my.list$ ggplot.string <- paste("ggplot(", quote(pdata), ",aes(sample =", match.term$WRES,
                                      ")) + stat_qq()", sep="")
      #cat(paste("code: ", my.list$ggplot.string, sep=""))
      generate.plot(lattice.call, my.list$lattice.string, ggplot.call, my.list$ggplot.string,
                    plot.dir, pdata, paste.string=FALSE)

      # Individual distribution of WRES:
      newstr <- paste("Individual distribiton of ", match.term$WRES, sep="")
      addCodeNote(newstr)

      ggplot.layout <- .pkplot$getHistGraph("others")$ind.layout      

      lattice.string <- paste("histogram(~", match.term$WRES, "|factor(", match.term$ID, "),",
                      "data=pdata, xlab=\"", match.term$WRES, "\")",
                     sep="")
      tmp.string1 <- paste("pdata[,\"", match.term$ID, "\"] <- factor(pdata[,\"", match.term$ID, "\"])", sep="")
      
      string.bin <- round(diff(range(pdata[,match.term$WRES]))/30,5)
      if (string.bin != 0)
      { 
          tmp.string2 <- paste("ggplot(pdata)+geom_histogram(aes(", match.term$WRES, 
                          "),", "binwidth=", string.bin,
                          ")",
                          "+facet_wrap(~", match.term$ID, ", ncol=", as.numeric(ggplot.layout[1]),
                          #"binwidth=", diff(range(pdata[,match.term$WRES]))/30, 
                          ")", sep="")
      }
      else
      {
          tmp.string2 <- paste("ggplot(pdata)+geom_histogram(aes(", match.term$WRES, 
                          "),", "binwidth=", diff(range(pdata[,match.term$WRES]))/30,
                          ")",
                          "+facet_wrap(~", match.term$ID, ", ncol=", as.numeric(ggplot.layout[1]),
                          #"binwidth=", diff(range(pdata[,match.term$WRES]))/30, 
                          ")", sep="")      
      }
      ggplot.string <- paste(tmp.string1, tmp.string2, sep=";")
      generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                    plot.dir, pdata, paste.string=FALSE)

      # Individual distribution of WRES (QQ):
      newstr <- paste("Individual distribiton of ", match.term$WRES, "(QQ)", sep="")
      addCodeNote(newstr)

      ggplot.layout <- .pkplot$getScatterGraph("others")$ind.layout

      lattice.string <- paste("qqmath(~", match.term$WRES, "|factor(", match.term$ID,
                     "), data=pdata",
                      ")",
                     #", prepanel = prepanel.qqmathline, panel = function(x, ...) {panel.qqmathline(x, ...);panel.qqmath(x, ...)})",
                     sep="")
      tmp.string1 <- paste("pdata[,\"", match.term$ID, "\"] <- factor(pdata[,\"", match.term$ID, "\"])", sep="")
      tmp.string2 <- paste("ggplot(", quote(pdata), ") + geom_point(aes(sample =",
                              match.term$WRES, "), stat = \"qq\")",
                              "+facet_wrap(~", match.term$ID, ", ncol=", as.numeric(ggplot.layout[1]),")",
                              sep="")
      ggplot.string <- paste(tmp.string1, tmp.string2, sep=";")

      generate.plot(lattice.call, lattice.string, 
        ggplot.call, ggplot.string, plot.dir, pdata, paste.string=FALSE)
       
    }

    ## global setup for scatter plot
    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")

    lattice.call <- "xyplot"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)
    
    # |WRES| vs PRED:
    tmp.string1 <- paste("pdata <- data.frame(pdata, WRES.ABS=abs(pdata[[\"",
                    match.term$WRES, "\"]]))", sep="")
    tmp.string2 <- paste("xyplot(", match.term$PRED, "~", "WRES.ABS",
                         ",data=pdata, type=c(\"p\", \"smooth\"))",
                         sep="")
    tmp.string3 <- paste("qplot(", "WRES.ABS", ",", match.term$PRED,
                         ",data=pdata, geom=c(\"point\", \"smooth\"), se=FALSE)",
                         sep="")
    lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
    ggplot.string <- paste(tmp.string1, tmp.string3, sep=";")
    
    newstr <- paste("|WRES| vs ", match.term$PRED,sep="")
    addCodeNote(newstr)
    generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                  plot.dir, pdata, paste.string=FALSE)

    #|WRES| vs PRED|Covariates:
    for (i in 1: length(match.term$COV))
    {
        newstr <- paste("|WRES| vs ", match.term$PRED, " | ", match.term$COV[i], sep="")
        addCodeNote(newstr)
        ggplot.layout <- .pkplot$getScatterGraph("others")$ind.layout

            tmp.id <- unique(pdata[,match.term$COV[i]])            
            if (length(tmp.id)>25) 
            {
                tmp.index <- which(pdata[,match.term$COV[i]] %in% tmp.id[1:25])
                tmp.pdata <- pdata[tmp.index,] 
            }
            else
            {
                tmp.pdata <- pdata
            }
        
        tmp.string1 <- paste("pdata <- data.frame(pdata, WRES.ABS=abs(pdata[[\"",
                    match.term$WRES, "\"]]))", sep="")
        tmp.string2 <- paste("xyplot(", match.term$PRED, "~ WRES.ABS|factor(", match.term$COV[i],
                     "), layout=c(", ggplot.layout[1], ",", ggplot.layout[2], "), data=pdata, as.table=TRUE)",
                     sep="")
        tmp.string3 <- paste("pdata[,\"", match.term$COV[i], "\"] <- factor(pdata[,\"", match.term$COV[i], "\"])", sep="")
        tmp.string4 <- paste("ggplot(pdata)+geom_point(aes(WRES.ABS,", match.term$PRED, "))",
                      "+facet_wrap(~", match.term$COV[i], ", ncol=", as.numeric(ggplot.layout[1]),")",
                      sep="")
        lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
        ggplot.string <- paste(tmp.string1, tmp.string3, tmp.string4, sep=";")
        
        generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                      plot.dir, tmp.pdata, paste.string=FALSE)

    }
    
    #|IWRES| vs IPRED|Covariates:
    # Autocorrelation of WRES: Plot of WRESi against WRESi+1 (autocorr.wres).
    
    newstr <- "WRESi.1 vs WRESi+1"
    addCodeNote(newstr)
    
    tmp.string1 <- paste("tmp.data <- pdata[[\"", match.term$WRES, "\"]]", sep="")
    tmp.string2 <- "auto.df <- data.frame(WRESi.0= tmp.data[-length(tmp.data)], WRESi.1=tmp.data[-1])"
    tmp.string3 <- "xyplot(WRESi.0~WRESi.1, data=auto.df, type=c(\"p\", \"smooth\"), xlab=\"WRESi.1\", ylab=\"WRESi.0\")"
    tmp.string4 <- "qplot(WRESi.1, WRESi.0, data=auto.df, geom=c(\"point\", \"smooth\"), se=FALSE, xlab=\"WRESi.1\", ylab=\"WRESi.0\")"

    lattice.string <- paste(tmp.string1, tmp.string2, tmp.string3, sep=";")
    ggplot.string <- paste(tmp.string1, tmp.string2, tmp.string4, sep=";")
    
    generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                  plot.dir, pdata, paste.string=FALSE)

    #Covariates vs |residuals| (bw): box-and-whisker plots of covariates plotted against absolute values of weighted residuals (abs.wres.vs.cov.bw).
    ## global boxplot setup
    lattice.call <- "bwplot"
    ggplot.call <- "qplot"

    for (i in 1: length(match.term$COV))
    {
        newstr <- paste(match.term$COV[i], " vs |RES|", "_bw", sep="")
        addCodeNote(newstr)
        
        tmp.string1 <- paste("pdata <- data.frame(pdata, RES.ABS=abs(pdata[[\"",
                    match.term$RES, "\"]]))", sep="")
        tmp.string2 <- paste("bwplot(RES.ABS ~ factor(", match.term$COV[i],
                        "),data=pdata)",
                        sep="")
        tmp.string3 <- paste("pdata[,\"", match.term$COV[i], "\"] <- factor(pdata[,\"",
                       match.term$COV[i], "\"])", sep="")
        tmp.string4 <- paste("ggplot(pdata)+geom_boxplot(aes(",
                        match.term$COV[i], ",RES.ABS))",
                        sep="")

        lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
        ggplot.string <- paste(tmp.string1, tmp.string3, tmp.string4, sep=";")
        
        generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                  plot.dir, pdata, paste.string=FALSE)
    }
    
    invisible(NULL)

}

# Parameters diagnostics
PKreport.6 <- function(pdata=deparse(substitute(realdata)))
{

    if (missing(pdata)) stop("Please input data!")
    ## require terms:
    match.term <- .pkplot$getTerm()
    if(is.null(match.term$PARA))
        stop(paste("No parameter", " is defined in var.name!", sep=""))

    ## global config
    general.list <- .pkplot$getGlobalConfig()
    
    ## create folder
    plot.dir <- general.list$para.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup
    lattice.list <- .pkplot$getHistGraph("lattice")
    ggplot.list <- .pkplot$getHistGraph("ggplot")

    lattice.call <- "histogram"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)
    
    # Distribution of parameters:
    if (length(match.term$PARA)> 0)
    {
      for (i in 1: length(match.term$PARA))
      {   
          newstr <- paste("Distribution of ", match.term$PARA[i], sep="")
          addCodeNote(newstr)
          ggplot.list.tmp <- ggplot.list

          bin.tmp <- round(diff(range(pdata[,match.term$PARA[i]]))/30,5)
          if (bin.tmp != 0) ggplot.list.tmp$binwidth <- bin.tmp
          else ggplot.list.tmp$binwidth <- diff(range(pdata[,match.term$PARA[i]]))/30
          
          my.list <- addList1(lattice.list, ggplot.list.tmp, match.term$PARA[i])
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, plot.dir, pdata)

      }

    #Distribution of parameters (QQ):
      for (i in 1: length(match.term$PARA))
      {
          newstr <- paste("Distribution of ", match.term$PARA[i], "(QQ)", sep="")
          addCodeNote(newstr)
          my.list$lattice.string <- paste("qqmath(~", match.term$PARA[i], ",data=", quote(pdata), ")", sep="")
          my.list$ggplot.string <- paste("ggplot(", quote(pdata), ") + geom_point(aes(sample =",
                                    match.term$PARA[i], "), stat = \"qq\")", sep="")
          #cat(paste("code: ", my.list$lattice.string, "\n", sep=""))
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, 
          plot.dir, pdata, paste.string=FALSE)

      }
    ## global setup for scatter plot
    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")

    lattice.call <- "xyplot"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)

    if (length(.pkplot$getTerm()$PARA)>1)
    {
        #Parameter vs parameter:
        for (i in 1: (length(match.term$PARA)-1))
        {
            newstr <- paste(match.term$PARA[i+1], " vs ", match.term$PARA[i], sep="")
            addCodeNote(newstr)
            my.list <- addList2(lattice.list, ggplot.list, match.term$PARA[i], match.term$PARA[i+1])
            generate.plot(lattice.call, my.list$lattice.string, 
            ggplot.call, my.list$ggplot.string, plot.dir, pdata)
        }

    #Scatterplot matrix of parameters

          newstr <- "Scatterplot matrix of Parameters"
          addCodeNote(newstr)

          tmp.string1 <- paste("para.data <- pdata[,", deparse(match.term$PARA), "]", sep="")
          tmp.string2 <- "splom(para.data)"
          tmp.string3 <- "plotmatrix(para.data)"
          my.list$lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
          my.list$ggplot.string <- paste(tmp.string1, tmp.string3, sep=";")
          
          generate.plot(lattice.call, my.list$lattice.string, ggplot.call, my.list$ggplot.string,
                        plot.dir, pdata, paste.string=FALSE)
      }
    
   }
    
}

## Covariate model
PKreport.7 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")
    ## require terms:
    match.term <- .pkplot$getTerm()
    require.term <- c("PARA", "ETA", "WRES", "COV")
    PK.match <- match(require.term, names(match.term))
    if(length(PK.match[is.na(PK.match)]) > 0)
        stop(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined in var.name!", sep=""))

    ## global config
    general.list <- .pkplot$getGlobalConfig()
    
    ## create folder
    plot.dir <- general.list$cov.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup for scatter plot
    lattice.list <- .pkplot$getScatterGraph("lattice")
    ggplot.list <- .pkplot$getScatterGraph("ggplot")

    lattice.call <- "xyplot"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)

    # Scatterplot matrix of covariates
    if (length(match.term$COV)>1)
    {
          newstr <- "Scatterplot matrix of Covariates"
          addCodeNote(newstr)

          tmp.string1 <- paste("cov.data <- pdata[,", deparse(match.term$COV)," ]", sep="")
          tmp.string2 <- "splom(cov.data)"
          tmp.string3 <- "plotmatrix(cov.data)"
          lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
          ggplot.string <- paste(tmp.string1, tmp.string3, sep=";")

          generate.plot(lattice.call, lattice.string, ggplot.call, ggplot.string,
                        plot.dir, pdata, paste.string=FALSE)
    }

    # Parameters vs covariates
    if ( (length(match.term$COV)>0) && (length(match.term$PARA) > 0) )
    {
      for (i in 1: length(match.term$PARA))
      {
          for (j in 1: length(match.term$COV))
          {
              newstr <- paste(match.term$PARA[i], " vs ", match.term$COV[j], sep="")
              addCodeNote(newstr)
              my.list <- addList2(lattice.list, ggplot.list, match.term$COV[j], match.term$PARA[i])
              generate.plot(lattice.call, my.list$lattice.string,
                            ggplot.call, my.list$ggplot.string, plot.dir, pdata)
          }

      }
     }

    # ETAs vs covariates:
    if ( (length(match.term$ETA)>0) && (length(match.term$COV) > 0) )
    {
      for (i in 1: length(match.term$ETA))
      {
          for (j in 1: length(match.term$COV))
          {
              newstr <- paste(match.term$ETA[i], " vs ", match.term$COV[j], sep="")
              addCodeNote(newstr)
              my.list <- addList2(lattice.list, ggplot.list, match.term$COV[j], match.term$ETA[i])
              generate.plot(lattice.call, my.list$lattice.string,
                            ggplot.call, my.list$ggplot.string, plot.dir, pdata)
          }

      }
    }

    # WRES vs covariates:
    if ( (length(match.term$WRES)>0) && (length(match.term$COV) > 0) )
    {
      for (i in 1: length(match.term$COV[i]))
      {
          newstr <- paste(match.term$WRES, " vs ", match.term$COV[i], sep="")
          addCodeNote(newstr)
          my.list <- addList2(lattice.list, ggplot.list, match.term$COV[i], match.term$WRES)
          generate.plot(lattice.call, my.list$lattice.string,
                        ggplot.call, my.list$ggplot.string, plot.dir, pdata)

      }
    }

}

# Random effects
PKreport.8 <- function(pdata=deparse(substitute(realdata)))
{
    if (missing(pdata)) stop("Please input data!")
    
    ## require terms:
    match.term <- .pkplot$getTerm()
    require.term <- c("ETA")
    PK.match <- match(require.term, names(match.term))
    if(length(PK.match[is.na(PK.match)]) > 0)
        stop(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined in var.name!", sep=""))

    ## global config
    general.list <- .pkplot$getGlobalConfig()

    ## create folder
    plot.dir <- general.list$eta.dir
    if(file.exists(plot.dir)) unlink(plot.dir, recursive=TRUE)
    if(!is.null(plot.dir)) dir.create(plot.dir)
    tmp.dir <- getwd()
    on.exit(setwd(tmp.dir))
    setwd(paste(tmp.dir, plot.dir, sep="/"))

    ## global setup for hist
    lattice.list <- .pkplot$getHistGraph("lattice")
    ggplot.list <- .pkplot$getHistGraph("ggplot")

    lattice.call <- "histogram"
    ggplot.call <- "qplot"

    lattice.list$data <- quote(pdata)
    ggplot.list$data <- quote(pdata)
    

    if ( (length(.pkplot$getTerm()$ETA)!=0) )
    {
      #Distribution of ETAS:
      for (i in 1: length(match.term$ETA))
      {
          newstr <- paste("Distribution of ", match.term$ETA[i], sep="")
          addCodeNote(newstr)
          ggplot.list.tmp <- ggplot.list
          
          bin.tmp <- round(diff(range(pdata[,match.term$ETA[i]]))/30, 5)
          if (bin.tmp != 0) ggplot.list.tmp$binwidth <- bin.tmp  
          else ggplot.list.tmp$binwidth <- diff(range(pdata[,match.term$ETA[i]]))/30
                  
          my.list <- addList1(lattice.list, ggplot.list.tmp, match.term$ETA[i])
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, plot.dir, pdata)

      }

      #Distribution of ETAs (QQ):
      for (i in 1: length(match.term$ETA))
      {
          newstr <- paste("Distribution of ", match.term$ETA[i], "(QQ)", sep="")
          addCodeNote(newstr)
          my.list$lattice.string <- paste("qqmath(~", match.term$ETA[i], ",data=", quote(pdata), ")", sep="")
          my.list$ggplot.string <- paste("ggplot(", quote(pdata), ") + geom_point(aes(sample =",
                                    match.term$ETA[i], "), stat = \"qq\")", sep="")
          generate.plot(lattice.call, my.list$lattice.string, 
          ggplot.call, my.list$ggplot.string, plot.dir, pdata, paste.string=FALSE)

      }
    
      #Scatterplot matrix of ETAs:
      if (length(.pkplot$getTerm()$ETA)>1)
      {
          newstr <- "Scatterplot matrix of ETA"
          addCodeNote(newstr)

          tmp.string1 <- paste("eta.data <- pdata[,", deparse(match.term$ETA), "]", sep="")
          tmp.string2 <- "splom(eta.data)"
          tmp.string3 <- "plotmatrix(eta.data)"
          my.list$lattice.string <- paste(tmp.string1, tmp.string2, sep=";")
          my.list$ggplot.string <- paste(tmp.string1, tmp.string3, sep=";")

          generate.plot(lattice.call, my.list$lattice.string, ggplot.call, my.list$ggplot.string,
                        plot.dir, pdata, paste.string=FALSE)
      }

    }
    

}


## Run diagnostics
#   input - pdata: data
#         - methods: 1-8: description is on the top
PKfigure <- function(pdata=deparse(substitute(realdata)), methods=c(1,2,3,4,5,6,7,8), clean=TRUE)
{  
    old.warn <- getOption("warn")
    options(warn=-1)
    on.exit(options(warn=old.warn))
          
    if (length(pdata)==0 || nrow(pdata)==0) stop("nonmemObj@tabdata is Not available!")
    
   if (length(methods)==0) stop("Please specify the methods you want to diagnose model!")

   if (max(methods) > 8 || min(methods) < 1) stop("There are only 8 methods ranging from 1 to 8!")
   if (clean) 
   {
      while (dev.cur()!=1) dev.off()
      PKclean()
   }

   general.list <- .pkplot$getGlobalConfig()
   if (1 %in% methods)
   {
      if(file.exists("univar")) unlink(general.list$univar, recursive=TRUE)
      if(file.exists("bivar")) unlink(general.list$bivar, recursive=TRUE)
      PKreport.1(pdata)
   }

   if (2 %in% methods)
   {
      if(file.exists("ind")) unlink(general.list$ind, recursive=TRUE)
      PKreport.2(pdata)
   }

   if (3 %in% methods)
   {
      if(file.exists("gof")) unlink(general.list$gof, recursive=TRUE)
      PKreport.3(pdata)
   }

   if (4 %in% methods)
   {
      if(file.exists("struct")) unlink(general.list$struct, recursive=TRUE)
      PKreport.4(pdata)
   }

   if (5 %in% methods)
   {
      if(file.exists("resid")) unlink(general.list$resid, recursive=TRUE)
      PKreport.5(pdata)
   }

   if (6 %in% methods)
   {
      if(file.exists("para")) unlink(general.list$para, recursive=TRUE)
      PKreport.6(pdata)
   }

   if (7 %in% methods)
   {
      if(file.exists("cov")) unlink(general.list$cov, recursive=TRUE)
      PKreport.7(pdata)
   }

   if (8 %in% methods)
   {
      if(file.exists("eta")) unlink(general.list$eta, recursive=TRUE)
      PKreport.8(pdata)
   }
   
   # set name for later to replace. Need to change.
   .pkplot$setDataName(deparse(substitute(pdata)))
   
}