glpenalty <-
function(X,x0 = NA,x1 = NA,p=NA,b=NA,type = NA,plotpenalty = TRUE,allowed.error = 0.005,invert = FALSE) {
  if (min(is.na(x0)) == FALSE & min(is.na(x1) == FALSE)) {
    glpenalty = vglpenalty(X,x0,x1,plotpenalty,allowed.error,invert)
  } else if (min(is.na(x1)) == FALSE & is.na(b) == FALSE & is.na(type) == FALSE) {
    glpenalty = bglpenalty(X,x1,b,type,plotpenalty,allowed.error,invert)
  } else if (is.na(p) == FALSE & is.na(b) == FALSE & is.na(type) == FALSE) {
    glpenalty = pglpenalty(X,p,b,type,plotpenalty,allowed.error,invert)
  } else {
    stop("arguments supplied did not match requirements")
  }
  glpenalty
}
