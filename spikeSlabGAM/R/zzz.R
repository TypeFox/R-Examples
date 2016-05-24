.onAttach <- function(libname, pkgname) {
  vers <- packageDescription("spikeSlabGAM")[["Version"]]
  packageStartupMessage("## ----- This is spikeSlabGAM ", vers, " ----- ##\n",
    "Please note that a recent update to gridExtra has made it necessary\n",
    "  to change the interface for <plot.spikeSlabGAM> starting in version 1.1-9.\n",
    "Instead of arguments 'rows', 'cols', 'widths', 'heights', 'maxPlotsPerPage',\n",
    "  it now accepts only 'nrow' and 'ncol'.\n",
    "Arguments 'widths' & 'heights' can still be defined and are handed over\n",
    "  to <gridExtra:::marrangeGrob>.\n",
    "Sorry for the inconvenience.")
}
