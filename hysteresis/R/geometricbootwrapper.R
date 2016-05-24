hystenv <- new.env()
hystenv$warningcountHysteresis <- 0
geometricbootwrapper <-
function(j,wr1,wr2,x.pred,y.pred,n,cbb,joint){
  tryCatch(geometricboot(j,wr1,wr2,x.pred,y.pred,n,cbb,joint),error=function (e) {
    hystenv$warningcountHysteresis <- hystenv$warningcountHysteresis + 1
    geometricbootwrapper(j,wr1,wr2,x.pred,y.pred,n,cbb,joint)
  })
      }
