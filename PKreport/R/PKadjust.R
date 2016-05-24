#############################################################################################
## File: PKadjust.R
## Author: Xiaoyong Sun
## Date: 10/12/2009
## Goal: Adjust PK plot
## Notes:
##      -
#############################################################################################

PKadjust <- function(figno, save=FALSE,...)
{
    adj.para <- list(...)
    pdata <- .pkplot$getPKData()
    
    mylist <- .pkplot$getPKCode(figno)
    mycall <- mylist[[length(mylist)]]
    mydir <- mylist[[length(mylist)-1]]
    mypara <- mylist[-c(length(mylist)-1, length(mylist))]

    if (length(adj.para)==0) print(do.call(mycall, mypara))
    else
    {
        for (i in 1:length(adj.para))
        {
            mypara[[names(adj.para)[i]]] <- adj.para[[i]]
        }
        print(do.call(mycall, mypara))
    }
    
    if (save)
    {
        myformat <- .pkplot$getConfig("save.format")
        save.call <- switch(myformat,
                            jpeg=jpeg,
                            bmp=bmp,
                            png=png,
                            tiff=tiff,
                            dev.new=dev.new)
        if(is.null(save.call)) save.call <- "pdf"
        
        tmp.dir <- getwd()
        on.exit(setwd(tmp.dir))
        setwd(paste(tmp.dir, mydir, sep="/"))
        
        fig.name <- addFigName(mydir, myformat)
        do.call(save.call, list(file=paste("adjust", fig.name, sep="_"), width=.pkplot$getConfig("width"),
                      height=.pkplot$getConfig("height")))
        print(do.call(mycall, mypara))
        dev.off()
    }

}

## TODO:
##    - reorder code as order: x, y, data=,...
PKcode <- function(filename="PKcode.txt")
{
    myfile <- filename
    PK.con <- file(myfile, "w")

    if (.pkplot$getPKCodeLen() == 0) stop("No data available to show! Please run PKfigure.")
    for (i in 1: .pkplot$getPKCodeLen())
    {
        one.code <- .pkplot$getPKCode(i)
        one.note <- .pkplot$getPKCodeNote(i)
        
        ## folder name
        writeLines(paste("\n#Folder: ", one.code[[2]], "   Figure ID: ", i, sep=""), con=PK.con, sep="\n")
        ## note
        writeLines(paste("#", one.note, sep=""), con=PK.con, sep="\n")
        ## code
        change.code <- gsub("pdata", .pkplot$getDataName(), one.code[[1]])
        writeLines(change.code, con=PK.con, sep="\n")
    }
    
    on.exit(close(PK.con))
}

PKnum <- function(exp.data)
{
    if( (!is.data.frame(exp.data)) && (!is.matrix(exp.data)) ) stop("The data must be data frame or matrix!")
    
    result <- matrix(NA, nrow=nrow(exp.data), ncol=ncol(exp.data))
    sapply(1:ncol(exp.data), function(i)
          {
              result[,i] <<- as.numeric(as.character(exp.data[,i]))
          })
    return(result)
}

