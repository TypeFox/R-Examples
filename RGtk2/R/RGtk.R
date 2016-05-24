.gtkInitCheck <- local({
  initialized <- NULL
  function(args="R")
  {
    if (is.null(initialized))
      initialized <<- .C("R_gtkInit", length(args), x=args, success = logical(1), 
        PACKAGE = "RGtk2")$success
    initialized
  }
})
gtkInit <- .gtkInitCheck

.gtkCleanup <- function() {
  .C("R_gtkCleanup", PACKAGE = "RGtk2")
}
