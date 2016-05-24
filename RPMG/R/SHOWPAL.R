`SHOWPAL` <-
  function(COLLIST, NAME=FALSE, NUM=FALSE, ncol=5, BACK="transparent")
  {
    if(missing(NAME)) { NAME=FALSE }
       if(missing(NUM)) { NUM=FALSE }
    if(missing(BACK)) { BACK=NULL }
    
    N  = length(COLLIST)
    
    if(missing(ncol))
      {

        if( (N%%2) == 0 )
          {
            ncol=6
          }
        else
          {
            ncol=5
          }
        
      }
    
    if( (N%%2) != (ncol%%2) )
      {
        ncol=ncol+1
      }
    
    ##  SHOWPAL(COLLIST)
    
    opar = par(no.readonly = TRUE)

    if(!is.null(BACK))
      {
        par(bg=BACK)
        bg = BACK
      }
    else
      {
        bg = opar$bg
      }
    
    par(mfrow=c(1,1))
    
    plot(c(0,1), c(0,1), type='n', axes=FALSE, xlab='', ylab='')
    
    nrow = ceiling(N/ncol)

    dx = 1/ncol
    dy =  1/nrow

###  print(c(ncol, nrow))
    if(NAME==FALSE)
      {
        LABS = rep("", times=length(COLLIST))
      }
    else
      {
        if(NUM)
          {
            LABS = paste(sep=":", seq(from=1, to=length(COLLIST)), COLLIST)
          }
        else
          {
            
            LABS = COLLIST
          }
      }
    
    YN = OPTREPLOT(LABS, ncol = ncol ,  cols =COLLIST, scol= bg, bcol= bg, cex=.8 , newplot=FALSE)

    title("Color Palette")
    par(opar)
    
    return(list(N=N, ncol=ncol, nrow=nrow, dx=dx, dy=dy))
  }

