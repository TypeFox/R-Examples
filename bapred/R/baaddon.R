baaddon <-
function(params, x, batch) {

  if(!(class(params) %in% c("fabatch", "combat", "svatrain", "meancenter",
    "standardize", "ratioa", "ratiog", "noba")))
    stop("Invalid class of 'params'.")
  
  if(class(params)=="fabatch") {
    return(fabatchaddon(params, x, batch))
  }
  if(class(params)=="combat") {
    return(combatbaaddon(params, x, batch))
  }
  if(class(params)=="svatrain") {
    return(svabaaddon(params, x))
  }
  if(class(params)=="meancenter") {
    return(meancenteraddon(params, x, batch))
  }
  if(class(params)=="standardize") {
    return(standardizeaddon(params, x, batch))
  }
  if(class(params)=="ratioa") {
    return(ratioaaddon(params, x, batch))
  }
  if(class(params)=="ratiog") {
    return(ratiogaddon(params, x, batch))
  }
  if(class(params)=="noba") {
    return(nobaaddon(params, x, batch))
  }

}
