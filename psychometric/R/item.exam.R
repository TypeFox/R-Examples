"item.exam" <- 
function (x, y = NULL, discrim = FALSE) 
{
    x <- na.exclude(as.matrix(x))
    if (!discrim) {
        discrim <- NA
        
    }
    else {
        discrim <- discrim(x)
          
    }
    
    k <- ncol(x)
    n <- nrow(x)
    TOT <- apply(x, 1, sum)
    TOT.woi <- TOT - (x)
    diff <- apply(x, 2, mean)
    rix <- cor(x, TOT, use = "complete")
    rix.woi <- diag(cor(x, TOT.woi, use = "complete"))
    sx <- apply(x, 2, sd)
    vx <- ((n - 1)/n) * sx^2
    if (is.null(y)) {
        riy <- NA
      }
    else {
        y <- y
        riy <- cor(x, y, use = "complete")
    }
    i.val <- riy * sqrt(vx)
    i.rel <- rix * sqrt(vx)
    i.rel.woi <- rix.woi * sqrt(vx)
    
    mat <- data.frame(Sample.SD = sx, Item.total = rix, Item.Tot.woi = rix.woi, Difficulty = diff, 
        Discrimination = discrim, Item.Criterion = riy, Item.Reliab = i.rel, Item.Rel.woi = i.rel.woi,
        Item.Validity = i.val)
    return(mat)
}
