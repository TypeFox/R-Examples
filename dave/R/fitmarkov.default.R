fitmarkov.default <-
function(veg,t,adjust=FALSE,...) {
     o.fitmarkov<- rfitmarkov(veg,t,adjust=FALSE)
     o.fitmarkov$call<- match.call()
     cat("Call:\n") 
     class(o.fitmarkov) <- "fitmarkov"
     print(o.fitmarkov$call)
     o.fitmarkov
  }
