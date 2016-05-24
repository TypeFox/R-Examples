
nlsbootwrapper <-
function(j,wr1,wr2,x.pred,y.pred,n,cbb,joint){
          tryCatch(nlsboot(j,wr1,wr2,x.pred,y.pred,n,cbb,joint),error=function (e) {
           hystenv$warningcountHysteresis <- hystenv$warningcountHysteresis + 1
            nlsbootwrapper(j=NULL,wr1,wr2,x.pred,y.pred,n,cbb,joint)
           })
      }
