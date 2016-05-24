## Open 'Life Expectancy' window

setLifeExp <-
function(){
  ## Check if window is active
  if (is.open("LE.win")){
    DALYfocus("LE.win")
  } else {
    setLifeExp.startup()
  }
}