fspa.default <-
function(veg,method,d.rev,n.groups,...) {
     o.fspa<- fspa2(veg,method,d.rev,n.groups)
     class(o.fspa)<- "fspa"
     o.fspa$call<- match.call()
     cat("Call:\n") 
     print(o.fspa$call)
     o.fspa
  }
