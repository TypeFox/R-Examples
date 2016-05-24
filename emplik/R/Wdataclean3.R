Wdataclean3 <- function(z, d, zc=rep(1,length(z)), wt = rep(1,length(z)) ) 
{  niceorder <- order(z,-d)
   sortedz <- z[niceorder]
   sortedd <- d[niceorder]
   sortedw <- wt[niceorder]
   sortedzc <- zc[niceorder]

   n <- length(sortedd)
   y1 <- sortedz[-1] != sortedz[-n]
   y2 <- sortedd[-1] != sortedd[-n]
   y3 <- sortedzc[-1] != sortedzc[-n]
   y <- y1 | y2 | y3

   ind <- c(which(y | is.na(y)), n)

   csumw <- cumsum(sortedw)

   list( value = sortedz[ind], dd = sortedd[ind],
         weight = diff(c(0, csumw[ind])) )
}
##############################################################
# this function sorts the data, and collapse them
# if there are true tie. and number of tie counted in weight.
# zc acts as a control of collapse: even if (z[i],d[i]) = (z[j],d[j])
# but zc[i] != zc[j] then obs. i and j will not collapse into one.
##############################################################
