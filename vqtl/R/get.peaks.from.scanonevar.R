#'  @title Get Local Maxima from Scanonevar
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @description \code{get.peaks.from.scanonevar} scans the genome for loci such that
#'    the locus to the left has a lower value and the locus to the right has a lower value.
#'    This value can be either LOD score or -log10(p-value)
#'
#'  @param x the \code{scanonevar} object from which the peaks will be identified
#'  @param thresh Optionally, the threshold over which a value has to be to be considered a peak.
#'    For example if one locus has a LOD score of 1 and the loci to its sides have LOD score of
#'    0.9, that's not really an interesting or "peaky" locus.  Defaults to 3 if x is in LOD units
#'    and 0.05 if x is in p-values
#'
#'  @return tbl_df of identified loci
#'
#'  @details none
#'

get.peaks.from.scanonevar <- function(x, thresh) {

  # hack to get R CMD CHECK to run without NOTEs that these globals are undefined
  full.lod <- mean.lod <- var.lod <- 'fake.global'
  emp.p.full.lod <- emp.p.mean.lod <- emp.p.var.lod <- 'fake.global'

  if (units(x) == 'lods') {

    if (missing(thresh)) { thresh <- 3 }

    peaks <- x %>%
      mutate(full.peak = (full.lod >= lag(full.lod, default = 0) & full.lod >= lead(full.lod, default = 0))) %>%
      mutate(mean.peak = (mean.lod >= lag(mean.lod, default = 0) & mean.lod >= lead(mean.lod, default = 0))) %>%
      mutate(var.peak = (var.lod >= lag(var.lod, default = 0) & var.lod >= lead(var.lod, default = 0)))
  }

  if (units(x) == 'emp.ps') {

    if (missing(thresh)) { thresh <- 0.05 }

    peaks <- x %>%
      mutate(full.peak = (emp.p.full.lod <= lag(emp.p.full.lod, default = 1) & emp.p.full.lod <= lead(emp.p.full.lod, default = 1))) %>%
      mutate(mean.peak = (emp.p.mean.lod <= lag(emp.p.mean.lod, default = 1) & emp.p.mean.lod <= lead(emp.p.mean.lod, default = 1))) %>%
      mutate(var.peak = (emp.p.var.lod <= lag(emp.p.var.lod, default = 1) & emp.p.var.lod <= lead(emp.p.var.lod, default = 1)))
  }

  return(peaks)
}