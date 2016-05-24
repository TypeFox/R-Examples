"ltraj2sldf" <-  function(ltr, byid = FALSE)
{
    ## Verifications
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class \"ltraj\"")

    ## Conversion
    lixy <- lapply(ltr,
                   function(x) Line(as.matrix(x[!is.na(x$x),c("x","y")])))
    id <- unlist(lapply(ltr, function(x) attr(x, "id")))
    bu <- unlist(lapply(ltr, function(x) attr(x, "burst")))

    if (byid) {
        lev <- levels(factor(id))
        re1 <- lapply(lev, function(x) Lines(lixy[id==x], ID=x))
        res <- SpatialLines(re1)
        df <- data.frame(id=lev)
        row.names(df) <- lev
    } else {
        res <- lapply(1:length(lixy),
                      function(i) Lines(list(lixy[[i]]), ID=bu[i]))
        res <- SpatialLines(res)
        df <- data.frame(id=id, burst=bu)
        row.names(df) <- bu
    }

    ## Output
    res <- SpatialLinesDataFrame(res, data=df)
    return(res)
  }

