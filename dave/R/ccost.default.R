ccost.default <-
function(veg,oldgr,newgr,y,...) {
     o.ccost<- ccost2(veg,oldgr,newgr,y=1)
     class(o.ccost)<- "ccost"
     o.ccost$call<- match.call()
     cat("Call:\n") 
     print(o.ccost$call)
     cat("cf=",o.ccost$ccost,"\n")
     o.ccost
 }
