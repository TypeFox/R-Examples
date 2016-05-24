"isfuncCOP" <-
function(cop=NULL, para=NULL, delta=0.002, ...) {

  # Therefore properties 2.2.2a-c of Nelsen (2006, p. 8) are investigated on a grid
  # where delta denotes the interval length. The code is broken in to five distinct
  # units so as the message() output can be used to inform the user. The last unit
  # is the 2-increasing and is very very likely the slowest "test".

  if(is.null(cop)) {
     warning("must have copula argument specified, returning NULL")
     return(NULL)
  }
  if(delta < 0 | delta > 0.5) {
     warning("invalid delta argument specified, returning NULL")
     return(NULL)
  }
  uvs <- seq(delta, 1-delta, by=delta)

  message("Checking condition C(u,0) = 0 for all u---", appendLF=FALSE)
  cond <- unique(cop(u=uvs, v=0, cop=cop, para=para, ...))
  if(length(cond) != 1 | cond[1] != 0) {
    message("FALSE.")
    return(FALSE)
  } else {
    message("TRUE.")
  }

  message("Checking condition C(0,v) = 0 for all v---", appendLF=FALSE)
  cond <- unique(cop(u=0, v=uvs, cop=cop, para=para, ...))
  if(length(cond) != 1 | cond[1] != 0) {
    message("FALSE.")
    return(FALSE)
  } else {
    message("TRUE.")
  }

  message("Checking condition C(u,1) = u for all u---", appendLF=FALSE)
  cond <- unique(uvs - cop(u=uvs, v=1, cop=cop, para=para, ...))
  if(length(cond) != 1 | cond[1] != 0) {
    message("FALSE.")
    return(FALSE)
  } else {
    message("TRUE.")
  }
  
  message("Checking condition C(1,v) = v for all v---", appendLF=FALSE)
  cond <- unique(uvs - cop(u=1, v=uvs, cop=cop, para=para, ...))
  if(length(cond) != 1 | cond[1] != 0) {
    message("FALSE")
    return(FALSE)
  } else {
    message("TRUE.")
  }

  message("Checking 2-increasing condition for all (u,v)---", appendLF=FALSE)
  tmp <- sapply(uvs, function(u) {
            if(any(densityCOP(u, uvs, cop=cop, para=para,
                              truncate.at.zero=FALSE, ...) < 0)) {
               message("FALSE.")
               return(FALSE)
            }})
  message("TRUE.")
  return(TRUE) # Copula properties hold in investigated points
}
