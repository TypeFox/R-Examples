recode <-
function(vec,original,modified) {
    if (length(original) != length(modified)) stop("original & modofied lengths differ")
    rvec<-vector(mode=mode(vec),length=length(vec))
    is.na(rvec)<-TRUE
    for (k in 1:length(original)) {
      rvec[vec==original[k]]<-modified[k]
    }
    return(rvec)
  }
