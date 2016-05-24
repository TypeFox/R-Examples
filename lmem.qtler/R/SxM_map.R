#' Name of the file containing genotypic (marker scores) information.
#'
#' The data is the well-known Steptoe x Morex doubled haploid population
#' developed in the early 90s by the North American Barley Mapping Project.
#' The objective was to improve in the understanding of the genetic basis of
#' agronomic and malting quality traits in barley.
#' The population consists of 150 doubled haploids lines;
#' of which 148 have been genotyped by SNP markers
#' (we use here 794 SNP markers).
#' The population was extensively evaluated for several agronomic and
#' malting quality traits (Hayes et al. 1993) in many locations and years
#' (US and Canada). In this example we use information on yield and heading
#' date in one of those trials.
#'
#' @format A data frame 794 row (markers) and 3 column.Column 1 gives the marker names,
#' column 2 the chromosome on which the marker has been mapped, and column 3 indicates
#' the position of the marker within the chromosome.
#'
#' @source {Hayes et al. 1993}
#'
"SxM_map"
