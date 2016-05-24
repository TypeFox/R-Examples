# ##
# Anything needed for the functioning of the namespace should 
# be handled at load/unload times by the .onLoad and .onUnload 
# hooks. For example, DLLs can be loaded (unless done by a 
#                                         useDynLib directive in 
#                                         the NAMESPACE file) 
# and initialized in .onLoad and unloaded in .onUnload. 
# Use .onAttach only for actions that are needed only when the 
# package becomes visible to the user (for example a start-up message) 
# or need to be run after the package environment has been created.

.onUnload <- function(libpath)
{
  library.dynam.unload("psd", libpath)
}
##
# executed after .onLoad is executed, once the namespace is visible to user
.onAttach <- function(...) {
  ##
  ## add some info to the environment
  psd::psd_envAssign("init", "initialized upon attach")
  ##
  packageStartupMessage(
    sprintf("Loaded psd (%s) -- Adaptive multitaper spectrum estimation",
            utils::packageVersion("psd")))
  ##
  ## options which can be reset by the user
  options(psd.ops=list(
    tapmin = 1,
    tapcap = 1000,
    names = list(
      fft = "working_fft",
      fft.padded = "fft_even_demeaned_padded",
      last.taper = "last_taper_sequence",
      last.psdcore = "last_psdcore_psd",
      last.psdcore.extrap = "last_psdcore_psd_extrap",
      series.even = "ser_orig_even",
      var.even = "ser_even_var",
      n.even = "len_even",
      n.even.half = "len_even_half",
      series.orig = "ser_orig",
      n.orig = "len_orig"
      )
  ))
}
