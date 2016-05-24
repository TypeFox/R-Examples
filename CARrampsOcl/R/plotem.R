plotem <- function(range1,range2,objname, range3, nrow,numcols,col,
rev.inds, breaks, title, sub)
{
        inds1 <- range1
        inds2 <- range2
        if(rev.inds[1])
           inds1 <- rev(inds1)
        if(rev.inds[2])
           inds2 <- rev(inds2)
        image.plot(range1, range2, 
            matrix(objname$y[range3], nrow = nrow)[inds1,inds2], 
            col = col, 
           breaks = breaks, xlab = "", ylab = "", main = title[1],
           sub=sub)
        image.plot(range1, range2, 
            matrix(objname$phi$phimean[range3], nrow = nrow)[inds1,inds2], 
            col = col, xlab = "", ylab = "", 
            breaks = breaks, main = title[2], sub=sub)
}
