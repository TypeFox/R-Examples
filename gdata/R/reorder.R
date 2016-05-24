reorder.factor <- function(x,
                           X,
                           FUN,
                           ...,
                           order=is.ordered(x),
                           new.order,
                           sort=mixedsort)
{
    constructor <- if (order) ordered else factor

    if(!missing(X) || !missing(FUN))
        {
            if(missing(FUN)) FUN <- 'mean'

            ## I would prefer to call stats::reorder.default directly,
            ## but it exported from stats, so the relevant code is
            ## replicated here:
            ## -->
            scores <- tapply(X = X, INDEX = x, FUN = FUN, ...)
            levels <- names(base::sort(scores, na.last = TRUE))
            if(order)
                ans <- ordered(x, levels=levels)
            else
                ans <- factor(x, levels=levels)
            attr(ans, "scores") <- scores
            ## <--
            return(ans)
        }
    else if (!missing(new.order))
      {
        if (is.numeric(new.order))
          new.order <- levels(x)[new.order]
        else
          new.order <- new.order
      }
    else
      new.order <- sort(levels(x))

    constructor(x, levels=new.order)
}
