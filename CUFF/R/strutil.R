# -*- coding: utf-8 -*-
###########################################################
### Author : Charles-Edouard Giguere                    ###
### Function to deal with strings.                      ###
###########################################################


### Overloading of paste with sep = "" using operator %+%.
'%+%' <- function(x,y){
    paste(x, y, sep = "")
}

### Operator %n% repeats string x, n times.
'%n%' <- function(x, y){
    sx=as.character(x)
    sy=as.integer(y)
    if(length(sx)>0 & length(sy)>0){
        mapply(function(x,y)paste(rep(x, y), collapse = ""),
        {
            if(length(sy) > 1 & length(sx) %in% 1){
                rep(sx, length(sy))
            } else{
                sx
            }
        },
        {
            if(length(sx) > 1 & length(sy) %in% 1){
                rep(sy, length(sx))
            } else{
                sy
            }
        },
        USE.NAMES = FALSE)
    }
    else{
        character()
    }
}


### Conversion of a numeric into a character using a formating.
### This is an alternative to sprintf.
### If nch is Null the function use the right amount to print
### the numerical variable with the right number of digits.
### If nch is not enough for the precision a star is returned.
numtostr=function(x,nch=NULL,digits=4){
    sx <- formatC(x,format="f",digits=digits)
    sx.length <- nchar(sx)
    if(!is.null(nch)){
        sx[sx.length>nch] <- (" " %n% (nch-1)) %+% "*"
        sx[sx.length<nch] <-  (" " %n% (nch-sx.length[sx.length<nch])) %+%
            sx[sx.length<nch]
    }
    sx
}



### Function that print a character with a specified format and align.
### if length is not enough 3 stars are printed padded with space.
strl <- function(x,length=max(nchar(x)),align="right"){
    x.n <- nchar(x)
    if(sum(x.n > length) > 0){
        if(align=="left")
            x[x.n>length] <-   "***"  %+% (" " %n% (length-3))
        else
            x[x.n>length] <- (" " %n% (length-3)) %+% "***"
    }
    if(sum(x.n<length)>0){
        if(align=="left")
            x[x.n<length] <- x[x.n<length] %+% (" " %n% (length-x.n[x.n<length]))
        else
            x[x.n<length] <- (" " %n% (length-x.n[x.n<length])) %+% x[x.n<length]
    }
    x
}
