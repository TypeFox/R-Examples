## Open 'options' window

DALYoptions <-
function(){
  ## Check if window is active
  if (is.open("opt.win")){
    DALYfocus("opt.win")
  } else {
    DALYoptions.startup()
  }
}