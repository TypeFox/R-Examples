.onAttach <- function(...) {
  if (interactive()) {
    OpenRSurvey()
  } else {
    msg <- "The RSurvey GUI is launched only in interactive sessions"
    packageStartupMessage(msg)
    return()
  }
  
}

.onLoad <- function(...) {
  if (interactive()) 
    LoadPackages()
  else
    return()
}
