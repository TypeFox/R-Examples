plot2Q <- function (objname, numcols=64, col=rev(terrain.colors(numcols)),
rev.inds = c(FALSE,FALSE))
{
    if (!(length(objname$n == 2)))
        print("Not a 2-dimensional problem; cannot plot with plot2D.")
    else {
       # require(fields)
       l1 <- objname$n[1]
       l2 <- objname$n[2]
       myrange <- range( c(objname$y, objname$phi$phimean) )
       mybreaks <-  seq(myrange[1], myrange[2],  length = numcols + 1)
       indices <- list( 1:l1, 1:l2)
       for( i in 1:2)
            if(rev.inds[i])
                indices[[i]] <- rev(indices[[i]])
                 
        par(mfcol = c(1, 2))
        image.plot(1:l1, 1:l2, matrix(objname$y, nrow = l1)[indices[[1]],
           indices[[2]]],
       #     col = rev(terrain.colors(numcols)),
            col=col,
            breaks = mybreaks,
            xlab = "", ylab = "", main = "Raw data")
        image.plot(1:l1, 1:l2, matrix(objname$phi$phimean, nrow = l1)[indices[[1]],
           indices[[2]]],
       #     col = rev(terrain.colors(numcols)),
            col=col,
            breaks = mybreaks,
            xlab = "", ylab = "",
            main = "Estimated underlying truth")
    }
}
