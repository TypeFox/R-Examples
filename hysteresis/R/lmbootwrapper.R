hystenv <- new.env()
hystenv$warningcountHysteresis <- 0
lmbootwrapper <-
function(j,wr1,wr2,x.pred,y.pred,n,cbb,joint){
  tryCatch(lmboot(j,wr1,wr2,x.pred,y.pred,n,cbb,joint),error=function (e) {
    hystenv$warningcountHysteresis <- hystenv$warningcountHysteresis + 1
    lmbootwrapper(j,wr1,wr2,x.pred,y.pred,n,cbb,joint)
  })
      }
