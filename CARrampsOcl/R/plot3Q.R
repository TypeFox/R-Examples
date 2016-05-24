plot3Q <- function( objname, numcols=64, col=rev(terrain.colors(numcols))
, plotdims=c(1,2), 
rev.inds=c(FALSE,FALSE), blocks=NULL, animate=FALSE, intv=3, title=c("Raw data","Estimated underlying truth"), sub=NULL ) 
{
# 
#  objname:  output object from CARramps.fit
#  numcols:  number of color shades
#  rev.cols:  if TRUE, reverse order of terrain.colors; if FALSE, use
#          default order
#  plotdims:  which 2 dimensions should be shown on each 2D plot
#  rev.inds: vector of two logicals.  If rev.inds[i] is TRUE, reverse the
#           indices in plotdims[i] before plotting
#  blocks:  which levels of the other dimension to loop over
#          NULL means all; otherwise specify as vector of indices
#  animate: if TRUE, there will be an animated paulse of length intv seconds
#          between plot displays; if FALSE, user is prompted to press enter
#          between plots
# title:  character vector of length 2; title[1] is title for first plot,
#          and title[2] for second plot.  
#          titles
    errors <- 0
    if(!(length(title) ==2))
    if (!(length(objname$n) == 3)) 
        {
         print("Not a 3-dimensional problem; cannot plot with plot3D.")
         errors <- 1
        }

    if( !(identical(plotdims,c(1,2)) | identical(plotdims,c(1,3)) | 
        identical(plotdims,c(2,3) )))  
        {
         print("Invalid plotdims.")
         errors <- 1
        }
    if( !errors )
    {

       #require(fields)    
       l <- objname$n
       if( is.null(blocks) )
           blocks <- 1:l[ (1:3)[-plotdims] ]  # 1 through length of other dim

#       mybreaks <- quantile(c(objname$y, objname$phi$phimean), 
#            seq(0, 1, length = numcols + 1))
       myrange <- range( c(objname$y, objname$phi$phimean) )
       mybreaks <-  seq(myrange[1], myrange[2],  length = numcols + 1)


       rangea <- list()
       for(i in 1:2)
            rangea[[i]] <- 1:l[plotdims[i]]

        nrow <- l[plotdims[1]]
        blocksize <- l[plotdims[1]] * l[plotdims[2]]
        if( identical( plotdims, c(1,3) ))
             oneblock <- rep((1:l[1]),l[3]) +
                     rep(  rep( (0:(l[3]-1)),each=l[1])* l[1] * l[2])
          if( identical( plotdims, c(2,3) ))
             oneblock <- seq( 1,(l[1]* l[2] * l[3]), by = l[1])

          par(mfrow=c(1,2))
 #    X11(width=8, height=6)


          for(j in blocks)
           {
            if( identical( plotdims, c(1,2) ))
               {
                rangea[[3]] <- ( (j-1) * blocksize + 1): (j * blocksize)
               }
            else if( identical( plotdims, c(1,3) ))
               { 
                rangea[[3]] <-  oneblock + (j-1) * l[1]
               }
            else if( identical( plotdims, c(2,3) ))
               {
                rangea[[3]] <- oneblock + (j-1)
               }
            if( !(is.null(sub)))
                cursub <- sub[j]
            else
                cursub <- NULL
            plotem(rangea[[1]],rangea[[2]],objname, rangea[[3]], nrow,numcols, 
                col, rev.inds, mybreaks, title, cursub)
            if( !(j==blocks[length(blocks)]))
               {
                if (animate)   # change MKC
                   {  
                     dev.flush()
                     Sys.sleep(intv)  ## pause for intv
                   }
               else
                  {
                    print("Press enter for next plot.")
                    scan(what = character(), sep = "\n", strip.white = T)
                  }
              }
            }  # for j in blocks
        
    } # if ! errors
 }
