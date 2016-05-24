if (is.R()) { ## R

"export.eps" <-
  function(FileName.in, Name.in="GSD2", ...) {
    ## Name.in is ignored in R
    dev.copy2eps(file=FileName.in, ...)
  }

} else { ## S-Plus
  
  export.eps <- function(FileName.in, Name.in="GSD2", ...) {
    ##  ... is ignored in S-Plus
    export.graph(FileName=FileName.in, Name=Name.in, ExportType = "EPS")
  }
  
}
