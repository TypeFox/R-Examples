powerDich <- function() {
  pbat.power(mode="dichotomous")
  return(invisible())
}
powerCont <- function() {
  pbat.power(mode="continuous")
  return(invisible())
}

fbatiLaunch <- function() {
  fbati()
  return(invisible())
}
fbatjLaunch <- function() {
  fbatj()
  return(invisible())
}
fbatcLaunch <- function() {
  fbatc()
  return(invisible())
}
fbatcStrategyStepLaunch <- function() {
  fbatcStrategyStep()
  return(invisible())
}

fbatgeLaunch <- function() {
  fbatge()
  return(invisible())
}

launchpad <- function() {
  ##library( fbati )  ## loads pbatR, fgui

  fguiNewMenu( c("pbatR","PBAT"), pbat )
  fguiNewMenu( c("pbatR","SEPARATOR") )
  fguiNewMenu( c("pbatR","Dichotomous trait power"), powerDich )
  fguiNewMenu( c("pbatR","Continuous trait power"), powerCont )

  fguiNewMenu( c("fbati","FBAT-I"), fbatiLaunch )
  fguiNewMenu( c("fbati","FBAT-J"), fbatjLaunch )
  fguiNewMenu( c("fbati","SEPARATOR") )
  fguiNewMenu( c("fbati","FBAT-C"), fbatcLaunch )
  fguiNewMenu( c("fbati","FBAT-C Stepwise"), fbatcStrategyStepLaunch )
  fguiNewMenu( c("fbati","SEPARATOR") )
  fguiNewMenu( c("fbati","FBAT-GxE"), fbatgeLaunch )
}
