## Open 'data' window

setData <-
function(n){
  ## Check if window is active
  win <- paste("data", n, sep = "")
  if (is.open(win)){
    DALYfocus(win)
  } else {
    setData.startup(n)
  }
}