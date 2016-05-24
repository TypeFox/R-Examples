imap <-
function (longlat = list(world.h.land, world.h.lake, world.h.island, world.h.pond.in.island, world.h.borders), 
    longrange, latrange, zoom = TRUE, col= c("black", "blue", "forestgreen", "dodgerblue", "cyan"), fill = TRUE, 
    poly = c("grey40", "blue", "forestgreen", "dodgerblue", NA), lwd = 1, keep.attr = TRUE, add.all = FALSE, bg = "grey81", tol = 0.05, ...) 
{
    par(bg = bg)

  
    if (is.matrix(longlat) | is.data.frame(longlat)) {
        LongLat <- longlat
        longlat <- list()
        longlat[[1]] <- list(ll=as.matrix(LongLat), col = col[1], lwd = lwd[1], poly=poly[1])
        rm(LongLat)
    }   

    if (!is.list(longlat[[1]])) {
    
        N <- length(longlat)
        col <- rep(col, length = N)
        lwd <- rep(lwd, length = N)
        poly <- rep(poly, length = N)
      
        LongLat <- longlat
        longlat <- list()
        for(i in 1:N)
             longlat[[i]] <- list(ll = LongLat[[i]], col=col[i], lwd=lwd[i], poly = poly[i])
        rm(LongLat)

     }  else {

        if(!keep.attr) {

           N <- length(longlat)
           col <- rep(col, length = N)
           lwd <- rep(lwd, length = N)
           poly <- rep(poly, length = N)
      
           LongLat <- longlat
           longlat <- list()
           for(i in 1:N)
               longlat[[i]] <- list(ll = LongLat[[i]]$ll, col=col[i], lwd=lwd[i], poly = poly[i])
           rm(LongLat)
        }
    }
   
    if (missing(longrange) | missing(latrange)) {
        max.ll <- apply(data.frame(lapply(longlat, function(x) apply(x[['ll']], 
            2, max, na.rm = TRUE))), 1, max)
        min.ll <- apply(data.frame(lapply(longlat, function(x) apply(x[['ll']], 
            2, min, na.rm = TRUE))), 1, min)
    }
    if (missing(longrange)) 
        longrange <- c(min.ll[1], max.ll[1])
    if (missing(latrange)) 
        latrange <- c(min.ll[2], max.ll[2])

    ll.out <- list()
    
    for (i in 1:length(longlat)) {

        ll.out[[i]] <- list(ll=matrix(

             imap.ll(longlat[[i]][['ll']], longrange, latrange, add = ifelse(add.all, TRUE, ifelse(i == 1, 
                FALSE, TRUE)), zoom = FALSE, col = longlat[[i]][['col']], poly = ifelse(fill, longlat[[i]][['poly']], NA),
                lwd = longlat[[i]][['lwd']], ...), 
              
        ncol = 2), col=longlat[[i]][['col']], lwd=longlat[[i]][['lwd']], poly = longlat[[i]][['poly']])
   }   
  
    if (zoom) {
        if (is.list(c1 <- locator(1))) {
            points(c1, pch = 3, cex = 3, col = 2)
            points(c2 <- locator(1), pch = 3, cex = 3, col = 2)
            if (abs(c1$x - c2$x) < tol & abs(c1$y - c2$y) < tol) {
                points(c1, pch = 0, cex = 3, col = 3)
                ll.out <- imap(longlat, zoom = TRUE, fill = fill, bg = bg, keep.attr = keep.attr, tol = tol, ...)
            }
            else ll.out <- imap(longlat, c(c1$x, c2$x), c(c1$y, 
                c2$y), zoom = TRUE, fill = fill, bg = bg, keep.attr = keep.attr, tol = tol, ...)
        }
    }

    for (i in length(ll.out):1) {

        if(all(is.na(ll.out[[i]]$ll)))
            ll.out[[i]] <- NULL   
    }        

    invisible(ll.out)
}

