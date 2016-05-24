## Helper function for 'setData()'
## Create combobox

defineCombo <-
function(frame, var, list){
  w <- 17
  if (length(list) == 4) w <- 12
  combo <- ttkcombobox(frame, state = "readonly", values = list,
                       textvariable = var, width = w)
  return(combo)
}
