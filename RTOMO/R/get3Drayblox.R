`get3Drayblox` <-
function(XNOD, YNOD, ZNOD , xo, yo, ztop, slowness=NULL)
{
  if(missing(slowness)) slowness=NULL
  ex = NULL
  why =  NULL
  zee =  NULL
  are =  NULL

  nnod = length(ZNOD)
  MIDZNOD = ZNOD[1:(nnod-1)] + (ZNOD[2:nnod]-ZNOD[1:(nnod-1)])/2


  FZNOD = findInterval( MIDZNOD    , ztop)

  for(k in 1:(length(XNOD)-1)  )
    {

      seglen = sqrt( (ZNOD[k]-ZNOD[k+1])^2+(XNOD[k]-XNOD[k+1])^2+(YNOD[k]-YNOD[k+1])^2 )
      if(seglen<=0) next

      IYZ = get2Drayblox(XNOD[k], YNOD[k], XNOD[k+1], YNOD[k+1], xo, xo , NODES=FALSE, PLOT=FALSE)


      ex = c(ex, IYZ$ix)
      why = c(why, IYZ$iy)
      zee = c(zee, rep(FZNOD[k], length(IYZ$ix)))

      cosalph = sqrt((XNOD[k]-XNOD[k+1])^2+(YNOD[k]-YNOD[k+1])^2 )/seglen
      if(cosalph<=0) 
        { are1 = IYZ$lengs }
      else
        {
          are1 = IYZ$lengs/cosalph
        }
      are = c(are, are1)

    }
  if(is.null(slowness)) 
    {
      tt = NULL
    }
  else
    {
      tt = sum(slowness[zee]*are)
    }
  return( list(ix=ex, iy=why, iz=zee, r=are, tt=tt) ) 
}

