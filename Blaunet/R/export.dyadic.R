export.dyadic <-
function(blauObj){
  if (is.null(blauObj$dyadic)){
    print("Nothing to export.")
  }
  else{
    return(blauObj$dyadic)
  }
}
