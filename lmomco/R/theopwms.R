"theopwms" <- 
function(para,nmom=5,verbose=FALSE) {
  if(nmom < 1) {
    warning("Number of PWMs requested is less than 1")
    return()
  }

  B <- vector(mode="numeric",length=nmom)
  for(r in seq(0,nmom-1)) { # for each  order of moment
    XofF <- function(F) {
      par2qua(F,para,paracheck=FALSE)*F^r
    }
    # Perform the numerical integration
    int <- integrate(XofF,0,1)
    B[r+1] <- int$value
    if(verbose == TRUE) { # Handy messages
      cat(c("abs.error=",int$abs.error,
            "  subdivisions=",int$subdivisions,
            "  message=",int$message,"\n"),sep="")
    }
  }

  z <- list(betas=B, source="theopwms")
  return(z)
}
