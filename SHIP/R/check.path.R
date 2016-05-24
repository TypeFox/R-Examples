check.path <-
function(p1,p2) {
  if ( !all(is.na(p1)) & !all(is.na(p2)) ) {
   for (i in 1:length(p1)) { if (p1[i] %in% p2) { return(1) } }
  }  
  return(0)
}

