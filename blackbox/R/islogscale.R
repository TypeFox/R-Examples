islogscale <- function(string,
                       scale=blackbox.getOption("FONKgScale"),
                       extraScale =blackbox.getOption("extraScale")
                       ) {
  if (length(string)>1L) return(sapply(string,islogscale,scale=scale,extraScale=extraScale))
  if (string %in% names(scale)) ## works fully only is string is a single string
    return((tolower(scale[string])=="logscale"))
  # ELSE
  if (string %in% names(extraScale)) ##
    # useful only for logscale plot of a non-kriging variable
    return((tolower(extraScale[string])=="logscale"))
  # ELSE at this point string is in no scale
  if (string=="latt2Ns2" && "condS2" %in% names(scale)) ##
    # for kriging on condS2 and LRT on latt2Ns2 => use the scale of condS2 for latt2Ns2 unless explicit opposite choice
    return(islogscale("condS2", scale=scale))
  # ELSE
  return(FALSE) ## string 'ln(L)' can occur in normal use...
}

