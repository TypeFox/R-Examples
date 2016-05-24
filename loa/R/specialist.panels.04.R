#in development code
#[TBC - NUMBER] functions 

#panel.YXCompare

#NOTE: much borrowed from lattice 

#to do
#callWithThis document or drop
#very minor


##############################
##############################
##panel.compareZcases
##############################
##############################

panel.compareZcases <- function (x=x, y=y, z=NULL, ..., loa.settings = FALSE) 
{
    if (loa.settings) 
        return(list(group.args = c("pch"), zcase.args = c(""),
                    ignore=c("col"), 
            default.settings = list(key = FALSE, grid = TRUE, 
                                    reset.xylims = c("refit.xylims", "zlim.in.ylim"))))
    extra.args <- list(...)
    if ("groups" %in% names(extra.args)) {
        if ("group.args" %in% names(extra.args) && length(extra.args$group.args) > 
            0) {
            temp <- as.numeric(factor(extra.args$groups, levels = extra.args$group.ids))
            for (i in extra.args$group.args) {
                extra.args[[i]] <- extra.args[[i]][temp]
            }
        }
        extra.args$groups <- NULL
    }
    if ("zcases" %in% names(extra.args)) {
        if ("zcase.args" %in% names(extra.args) && length(extra.args$zcase.args) > 
            0) {
            temp <- as.numeric(factor(extra.args$zcases, levels = extra.args$zcase.ids))
            for (i in extra.args$zcase.args) {
                extra.args[[i]] <- extra.args[[i]][temp]
            }
        }
########extra.args$zcases <- NULL
    }
    if (isGood4LOA(extra.args$grid)) 
        panel.loaGrid(panel.scales = extra.args$panel.scales, 
            grid = extra.args$grid)
    extra.args$grid <- NULL

    if(!is.null(z)){
       if(!"zcases" %in% names(extra.args))
           extra.args$zcases <- rep("default", length(z))
       if(!"zcase.ids" %in% names(extra.args))
           extra.args$zcase.ids <- "default"
    }
    temp <- length(extra.args$zcase.ids) + 1

#note reversal of ids
#might want to make this a global change?

    if("zcase.ids" %in% names(extra.args))
         extra.args$zcase.ids <- rev(extra.args$zcase.ids)
    
    extra.args$col <-if("col" %in% names(extra.args))
                     rep(extra.args$col, length.out=3) else 
                     rev(do.call(colHandler, listUpdate(extra.args, list(z=1:(temp+1), ref=1:(temp+1), zlim=NULL)))[-1])
    if("line.col" %in% names(extra.args))
           extra.args$col[1] <- extra.args$line.col

    if(temp>1){
        x <- x[extra.args$zcases == extra.args$zcase.ids[1]]
        y <- y[extra.args$zcases == extra.args$zcase.ids[1]]
        x1 <- c(x, rev(x))
        for(i in (temp-1):1){
            y1 <- z[extra.args$zcases == extra.args$zcase.ids[i]]
            y1 <- if(i>1)
                    c(y1, rev(z[extra.args$zcases == extra.args$zcase.ids[i-1]]))
                       else c(y1, rev(y))
            do.call(panel.polygon, listUpdate(extra.args, list(x=x1, y=y1, col=extra.args$col[i+1], alpha=0.1, border=FALSE)))
        }
        for(i in (temp-1):1){
            y1 <- z[extra.args$zcases == extra.args$zcase.ids[i]]
            do.call(panel.lines, listUpdate(extra.args, list(x=x, y=y1, col=extra.args$col[i+1])))
        }
    }
    do.call(panel.lines, listUpdate(extra.args, list(x=x, y=y, col=extra.args$col[1])))
}

