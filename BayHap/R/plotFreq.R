`plotFreq` <-
function(res,acf=NULL,d=NULL,rmean=NULL,tr=NULL){

          keep.rares.plot<-FALSE

          if (is.null(c(acf,d,rmean,tr))) {
             plotACF(res,keep.rares.acf=keep.rares.plot)
             plotDensity(res,keep.rares.density=keep.rares.plot)
             plotRmean(res,keep.rares.rmean=keep.rares.plot)
             plotTrace(res,keep.rares.trace=keep.rares.plot)
          }

          if ((!is.null(acf))&&(acf)) plotACF(res,keep.rares.acf=keep.rares.plot)
          if ((!is.null(d))&&(d))     plotDensity(res,keep.rares.density=keep.rares.plot)
          if ((!is.null(rmean))&&(rmean))   plotRmean(res,keep.rares.rmean=keep.rares.plot)
          if ((!is.null(tr))&&(tr))   plotTrace(res,keep.rares.trace=keep.rares.plot)


}

