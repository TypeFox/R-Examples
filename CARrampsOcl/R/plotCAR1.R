plotCAR1 <- function( objname, numcols=64, col=rev(terrain.colors(numcols)),
cardims,
rev.inds=c(FALSE,FALSE), title=c("Raw data","Estimated underlying truth"), sub=NULL ) 
{
# 
#  objname:  output object from CARramps.fit
#  numcols:  number of color shades
#  col:  color palette for image; if default is overridden, then 
#        numcols must have same explicit numerical value as in numcols argument
#  plotdim:   CAR1 matrix in Q list
# cardims:  number of rows and colums of data in CAR1; same as in "content"
#          part of Q list
#  rev.inds: vector of two logicals.  If rev.inds[i] is TRUE, reverse the
#           indices in plotdim[i] before plotting
#  blocks:  which levels of the other dimension to loop over
#          NULL means all; otherwise specify as vector of indices
#  animate: if TRUE, there will be an animated paulse of length intv seconds
#          between plot displays; if FALSE, user is prompted to press enter
#          between plots
# title:  character vector of length 2; title[1] is title for first plot,
#          and title[2] for second plot.  
#          titles
    
    plotdim <- 1
    blocks <- 1
    errors <- 0
    if(!(length(title) ==2))
    if ((length(objname$n) == 1)) 
        {
         print("Not a CAR1; cannot plot with plotCAR1.")
         errors <- 1
        }

    if( !errors )
    {
       #require(fields)
       l <- objname$n
       if( is.null(blocks) )
           blocks <- 1:l[ (1:2)[-plotdim] ]  # 1 through length of other dim

       #mybreaks <- quantile(c(objname$y, objname$phi$phimean), 
       #     seq(0, 1, length = numcols + 1))
       myrange <- range( c(objname$y, objname$phi$phimean) )
       mybreaks <-  seq(myrange[1], myrange[2],  length = numcols + 1)

       rangea <- list()
       for(i in 1:2)
            rangea[[i]] <- 1:cardims[i]

        nrow <- cardims[1]
        blocksize <- cardims[1] * cardims[2]

          par(mfrow=c(1,2))
#     X11(width=8, height=6)


          for(j in blocks)
           {
            rangea[[3]] <- ( (j-1) * blocksize + 1): (j * blocksize)
            if( !(is.null(sub)))
                cursub <- sub[j]
            else
                cursub <- NULL
            plotem(rangea[[1]],rangea[[2]],objname, rangea[[3]], nrow,numcols, 
                col, rev.inds, mybreaks, title, cursub)
            }  # for j in blocks
        
    } # if ! errors
 }
