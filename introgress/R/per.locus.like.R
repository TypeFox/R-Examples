per.locus.like <-
function(a1, a2, locustype, r1, r2, s1, s2, h){
  p <- h * r1 + (1 - h) * s1
  q <- h * r2 + (1 - h) * s2
  if(locustype== "C" | locustype=="c"){
    if (length(r1)==0 | length(r2)==0 | length(s1)==0 | length(s2)==0){
      return(NA)
    }
    else if(r1 == s1 && r2 == s2){
      return(NA)
    }
  }
  else if (locustype=="H" | locustype=="h" | locustype=="D" | locustype=="d"){
    if (length(r1)==0 | length(s1)==0) return(NA)
    else if (r1 == s1) return(NA)
  }
  if (locustype == "D" | locustype=="d"){
    if (is.na(a1)==FALSE){
      if(a1 == 1) { ## band present
        return(log(p^2 + 2 * p * q))
      }
      else if (a1 == 0){ ## band absent
        return(log(p^2))  ## r1 and s1 will be for the state that is observed
      }
    }
    else return(NA)
  }
  else if (locustype == "H" | locustype== "h"){  ## haploid
    return(log(p))
  }
  ## if data are codominant
  else if (is.na(a1)==FALSE & is.na(a2)==FALSE){
    if(a1 == a2){
      return(log(p^2))
    }
    else if (a1 != a2){
      return(log(2 * p * q))
    }
  }
  else {
    return(NA)
  }
}

