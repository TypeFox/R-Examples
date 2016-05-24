# split.labs.R: functions for generating split.labels

split.labs.wrapper <- function(x, split.fun, split.fun.name,
                               split.prefix, split.suffix,
                               right.split.prefix, right.split.suffix,
                               type, clip.left.labs, clip.right.labs, xflip, digits,
                               varlen, faclen, facsep, eq, lt, ge)
{
    labs <- internal.split.labs(x, type, clip.left.labs, clip.right.labs, xflip,
                                digits, varlen, faclen, facsep, eq, lt, ge,
                                split.prefix, right.split.prefix,
                                split.suffix, right.split.suffix)

    if(!identical(split.fun.name, "internal.split.labs")) { # call user's split.fun?
        check.func.args(split.fun, "split.fun",
                        function(x, labs, digits, varlen, faclen) NA)
        labs <- split.fun(x, labs, digits, varlen, faclen)
    }
    # check returned labs because split.fun may be user supplied
    check.returned.labs(x, labs, split.fun.name)
    labs
}
# Modified version of labels.rpart.
# This uses format0 instead of formatg and has various other extensions.
internal.split.labs <- function(x, type, clip.left.labs, clip.right.labs, xflip,
                                digits, varlen, faclen, facsep, eq, lt, ge,
                                split.prefix, right.split.prefix,
                                split.suffix, right.split.suffix)
{
    frame <- x$frame
    if(nrow(frame) == 1) # special case, no splits?
        return("root")   # NOTE: return
    is.leaf <- is.leaf(frame)
    split.var.names <- frame$var[!is.leaf]  # variable names for the (primary) splits
    split.var.names <- as.character(split.var.names) # factor levels to character
    clip.left.labs  <- recycle(clip.left.labs,  split.var.names)
    clip.right.labs <- recycle(clip.right.labs, split.var.names)

    # isplit is the row index of the primary split in x$splits
    index <- cumsum(c(1, frame$ncompete + frame$nsurrogate + !is.leaf))
    isplit  <- index[c(!is.leaf, FALSE)]

    split <- get.lsplit.rsplit(x, isplit, split.var.names,
                               type, clip.left.labs, clip.right.labs, xflip,
                               digits, faclen, facsep, eq, lt, ge)

    # We now have something like this:
    #    split.var.names:   sex   age     pclass     sibsp   pclass
    #    split$lsplit:      mal   >=9.5   =2nd,3rd   >=2.5   =3rd
    #    split$rsplit:      fml   <9.5    =1st       <2.5    =1st,2nd

    paste.split.labs(frame, split.var.names, split$lsplit, split$rsplit,
                     type, clip.left.labs, clip.right.labs, xflip, varlen,
                     split.prefix, right.split.prefix,
                     split.suffix, right.split.suffix)
}
get.lsplit.rsplit <- function(x, isplit, split.var.names,
                              type, clip.left.labs, clip.right.labs, xflip,
                              digits, faclen, facsep, eq, lt, ge)
{
    frame <- x$frame
    is.leaf <- is.leaf(frame)
    ncat  <- x$splits[isplit, "ncat"]
    lsplit <- rsplit <- character(length=length(isplit))
    is.con <- ncat <= 1             # continuous vars (a logical vector)
    if(any(is.con)) {               # any continuous vars?
        cutpoint <- format0(x$splits[isplit[is.con], "index"], digits)
        is.less.than <- ncat < 0
        lsplit[is.con] <- paste0(ifelse(is.less.than, lt, ge)[is.con], cutpoint)
        rsplit[is.con] <- paste0(ifelse(is.less.than, ge, lt)[is.con], cutpoint)
    }
    is.cat <- ncat > 1              # categorical variables (a logical vector)
    if(any(is.cat)) {               # any categorical variables?
        # jrow is the row numbers of factors within lsplit and rsplit
        # cindex is the index on the "xlevels" list
        jrow <- seq_along(ncat)[is.cat]
        crow <- x$splits[isplit[is.cat], "index"] # row number in csplit
        xlevels <- attr(x, "xlevels")
        cindex <- match(split.var.names, names(xlevels))[is.cat]
        # decide if we must add a "=" prefix
        paste.left.eq  <- !is.fancy(type) | (if(xflip) !clip.right.labs else !clip.left.labs)
        paste.right.eq <- !is.fancy(type) | (if(xflip) !clip.left.labs  else !clip.right.labs)
        for(i in 1:length(jrow)) {
            node.xlevels <- my.abbreviate(xlevels[[cindex[i]]],
                                          faclen, one.is.special=TRUE)
            j <- jrow[i]
            splits <- x$csplit[crow[i],]
            # splits is 1=left 2=neither 3=right
            left  <- (1:length(splits))[splits == 1]
            right <- (1:length(splits))[splits == 3]
            collapse <- if(faclen==1) "" else facsep
            lsplit[j] <- paste(node.xlevels[left],  collapse=collapse)
            if(paste.left.eq[i])
                lsplit[j] <- paste0(eq, lsplit[j])
            rsplit[j] <- paste0(node.xlevels[right], collapse=collapse)
            if(paste.right.eq[i])
                rsplit[j] <- paste0(eq, rsplit[j])
        }
    }
    list(lsplit=lsplit, rsplit=rsplit)
}
# Paste the various constituents to create the split labels vector.
# On entry we have something like this:
#    split.var.names sex   age     pclass     sibsp   pclass
#    lsplit:         mal   >=9.5   =2nd,3rd   >=2.5   =3rd
#    rsplit:         fml   <9.5    =1st       <2.5    =1st,2nd

paste.split.labs <- function(frame, split.var.names, lsplit, rsplit,
                             type, clip.left.labs, clip.right.labs, xflip, varlen,
                             split.prefix, right.split.prefix,
                             split.suffix, right.split.suffix)
{
    nodes <- as.numeric(row.names(frame))
    is.right <- nodes %% 2 == 1
    is.leaf <- is.leaf(frame)
    parent <- match(nodes %/% 2, nodes[!is.leaf])
    split.var.names <- my.abbreviate(split.var.names, varlen)
    left.names <- right.names <- split.var.names
    if(is.fancy(type)) {
        if(xflip)
            right.names[clip.left.labs] <- ""
        else
            left.names[clip.left.labs] <- ""
        if(xflip)
            left.names[clip.right.labs] <- ""
        else
            right.names[clip.right.labs] <- ""
    }
    if(is.null(right.split.prefix))
        right.split.prefix <- split.prefix
    if(is.null(right.split.suffix))
        right.split.suffix <- split.suffix

    # the heart of this function

    labs  <- paste0(split.prefix, left.names[parent], lsplit[parent], split.suffix)
    labs[1] <- "root" # was NA NA
    labs[is.right] <- paste0(right.split.prefix,
                             right.names[parent],
                             rsplit[parent],
                             right.split.suffix)[is.right]
    labs
}
