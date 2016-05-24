#' Pause in interactive mode and continue on <Enter>.
#' @export
pause = function() {
  if (interactive())
    readline("Pause. Press <Enter> to continue.")
}
