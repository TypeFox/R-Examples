`text2spc.fnc` <-
function(text) {
  if (!require("zipfR", quietly = TRUE)) {
    stop("please install the zipfR library first")
  } else {
    tab = table(table(text))
    return(spc(m = as.numeric(names(tab)), Vm = as.numeric(tab)))
  }
}

