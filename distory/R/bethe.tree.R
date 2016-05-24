bethe.tree <- function(tips, level.lengths = NULL, outgroup="O", outgroup.dist=1)
{
    if(length(tips) != length(unique(tips)))
    {
        stop("Not all tips are unique.\n")
    }

    if(!is.null(level.lengths))
    {
        if(length(level.lengths) == 1)
        {
            level.lengths = rep(level.lengths, 2)
        }
    }

    if(is.null(level.lengths) || length(level.lengths) == 0)
    {
        level.lengths = 1
    }

    groupings <- as.numeric(sapply(1:(length(tips)/2), function(x) rep(x,2)))

    d <- level.lengths[1]
    level.lengths <- level.lengths[-1]
    tips <- lapply(split(tips,groupings), function(x) { sprintf("(%s:%d,%s:%d)",x[1],d,x[2],d)} ) 
    if(length(tips) == 1)
    {
        read.tree(text=sprintf("(%s:%d, %s:%d);", tips[[1]], d, outgroup, outgroup.dist))
    }
    else
    {
        bethe.tree(tips, level.lengths, outgroup, outgroup.dist)
    }
}

