.onAttach <- function(library, pkg)
{
  # require("stats4")  
  # require("methods")
  # require("mnormt")
  # require("numDeriv")
  if(interactive())
  {
    # pkg <- Package("sn")
    meta <- packageDescription("sn")
    packageStartupMessage(
      "Package 'sn', ", meta$Version, " (", meta$Date, "). ",
      "Type 'help(SN)' for summary information.\n",
      "The package redefines function 'sd' but its usual working is unchanged.")
  }
  invisible()
}
