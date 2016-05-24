## Open 'population' window

setPop <-
function(){
  ## Check if window is active
  if (is.open("pop.win")){
    DALYfocus("pop.win")
  } else {
    setPop.startup()
  }
}