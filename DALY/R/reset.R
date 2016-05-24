## Close all active windows

reset <-
function(){
  ## Close any open windows
  if (DALYexists("pop.win")) tkdestroy(DALYget("pop.win"))
  if (DALYexists("LE.win")) tkdestroy(DALYget("LE.win"))
  for (i in seq(8))
    if (DALYexists(paste("data", i, sep = "")))
      tkdestroy(DALYget(paste("data", i, sep = "")))

  ## Reset all variables, *except* 'LE'
  resetVariables()
}