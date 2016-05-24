areDataNotLoaded <-
function() {
  tmp <- ls(pos = 1)
  if(sum(tmp == "geData") < 1 ) {
    message("The geData matrix is not found. Please Load!")
    return(TRUE)
  }
  if(sum(tmp == "stData") < 1 ) {
    message("The stData survival structure is not found. Please Load!")
    return(TRUE)
    }
  return(FALSE)
}
