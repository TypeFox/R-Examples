ba <-
function(x, y, batch, method=c("fabatch", "combat", "sva", "meancenter", "standardize", "ratioa", "ratiog", "none"), ...) {

  if(!(method %in% c("fabatch", "combat", "sva", "meancenter",
    "standardize", "ratioa", "ratiog", "none")))
    stop("Input parameter 'method' has to be one of the following:\n'fabatch',  'combat', 'fsva', 'meancenter', 'standardize', 'ratioa', 'ratiog', 'none'.")

  if(method=="fabatch") {
    return(fabatch(x, y, batch, ...))
  }
  if(method=="combat") {
    return(combatba(x, batch))
  }
  if(method=="sva") {
    return(svaba(x, y, batch, ...))
  }
  if(method=="meancenter") {
    return(meancenter(x, batch))
  }
  if(method=="standardize") {
    return(standardize(x, batch))
  }
  if(method=="ratioa") {
    return(ratioa(x, batch))
  }
  if(method=="ratiog") {
    return(ratiog(x, batch))
  }
  if(method=="none") {
    return(noba(x, batch))
  }  
  
}
