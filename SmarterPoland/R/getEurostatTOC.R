getEurostatTOC <-
function() {
  setEurostatTOC()
  invisible(get(".eurostatTOC", envir = .SmarterPolandEnv))
}
