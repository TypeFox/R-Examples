.onAttach <- function(library, pkg)
{
  # Rv <- R.Version()
  # if(Rv$major < 2 |(Rv$major == 2 && Rv$minor < 2.0))
  #  stop("This package requires R 2.2.0 or later")
  if(interactive())
  {
    meta <- packageDescription("circular")
    packageStartupMessage(
         "Package 'circular', ", meta$Version, " (", meta$Date, "). ",
         "Type 'help(Circular)' for summary information. \n Please report any bugs or comments to <Claudio Agostinelli> claudio@unive.it \n The package redefine how function 'var' and 'sd' works \n In particular, (try 'methods(var)' and 'methods(sd)') \n notice that 'var.default' ('sd.default')is an alias for the original 'var' ('sd') function \n and that a method for data.frame is available.\n")
  }
  invisible()
}
