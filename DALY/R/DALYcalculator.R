## Open DALY Calculator main window

DALYcalculator <-
function(){
  ## Check if window is active
  if (is.open("main")){
    DALYfocus("main")
  } else {
    DALYcalculator.startup()
  }
}