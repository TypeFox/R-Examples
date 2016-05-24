#' @title Process Control and Internal Validation of Forensic STR Kits
#' @docType package
#' @name strvalidator-package
#' @author Oskar Hansson \email{oskar.hansson@@fhi.no}
#' @keywords package
#' @description STR-validator is a free and open source R-package intended for
#' process control and internal validation of forensic STR DNA typing kit.
#' Its graphical user interface simplifies the analysis of data exported from
#' e.g. GeneMapper software, without extensive knowledge about R. It provides 
#' functions to import, view, edit, and export data. After analysis the
#' results, generated plots, heat-maps, and data can be saved in projects for
#' easy access. Currently, analysis modules for stutter, balance, dropout,
#' mixture, concordance, typing result, precision, pull-up, analytical
#' thresholds, and drop-in are available. In addition there are functions
#' to analyse the GeneMapper bins- and panels files. EPG like plots can be
#' generate from data. STR-validator can greatly increase the speed of
#' validation by reducing the time and effort needed to analyse the data.
#' It allows exploration of the characteristics of DNA typing kits according
#' to ENFSI and SWGDAM recommendations (see references). This facilitates the
#' implementation of probabilistic interpretation of DNA results.\cr\cr
#' 
#' Effort has been made to assure correct results. Refer to the main website
#' for a list of functions specifically tested at build time.\cr\cr
#' 
#' Click \code{Index} at the bottom of the page to see a complete list of
#' functions.\cr\cr
#' 
#' Created and maintained by:\cr
#' Oskar Hansson, Department of Forensic Biology (NIPH, Norway)\cr\cr
#' 
#' General user information and tutorials:\cr
#' \url{https://sites.google.com/site/forensicapps/strvalidator}\cr\cr
#' 
#' Facebook user community:\cr
#' \url{https://www.facebook.com/pages/STR-validator/240891279451450?ref=tn_tnmn}\cr\cr
#' \url{https://www.facebook.com/groups/strvalidator/}\cr\cr
#' 
#' Please report bugs to:\cr
#' \url{https://github.com/OskarHansson/strvalidator/issues}\cr\cr
#' 
#' The source code and collaborative community is hosted at GitHub:\cr
#' \url{https://github.com/OskarHansson/strvalidator}\cr\cr
#' 
#' @references
#' Recommended Minimum Criteria for the Validation of Various Aspects of the DNA Profiling Process
#' \url{http://www.enfsi.eu/sites/default/files/documents/minimum_validation_guidelines_in_dna_profiling_-_v2010_0.pdf}
#' Validation Guidelines for Forensic DNA Analysis Methods (2012)
#' \url{http://media.wix.com/ugd/4344b0_cbc27d16dcb64fd88cb36ab2a2a25e4c.pdf}
#' 
NULL

#' ESX17 Positive Control Profile
#' 
#' A dataset in 'GeneMaper' format containing the DNA profile of
#' the ESX17 positive control sample with homozygotes as one entry.
#' 
#' @docType data
#' @keywords datasets
#' @name ref1
#' @usage data(ref1)
#' @format A data frame with 17 rows and 4 variables
NULL

#' SGMPlus example data
#' 
#' A slimmed reference dataset containing an arbitrary SGMPlus DNA profile.
#' 
#' @docType data
#' @keywords datasets
#' @name ref2
#' @usage data(ref2)
#' @format A data frame with 16 rows and 3 variables
NULL

#' ESX17 Positive Control Profile
#' 
#' A dataset in 'GeneMaper' format containing the DNA profile of
#' the ESX17 positive control sample with homozygotes as two entries.
#' 
#' @docType data
#' @keywords datasets
#' @name ref11
#' @usage data(ref11)
#' @format A data frame with 17 rows and 4 variables
NULL


#' Typing data in 'GeneMapper' format
#' 
#' A dataset containing ESX17 genotyping result for 8 replicates
#' of the positive control sample, a negative control and ladder.
#' 
#' @docType data
#' @keywords datasets
#' @name set1
#' @usage data(set1)
#' @format A data frame with 170 rows and 13 variables
NULL

#' SGMPlus example data
#' 
#' A slimmed dataset containing SGM Plus genotyping result for 2 replicates
#' of 'sampleA'.
#' 
#' @docType data
#' @keywords datasets
#' @name set2
#' @usage data(set2)
#' @format A data frame with 32 rows and 5 variables
NULL

#' ESX17 example data for dropout analysis.
#' 
#' Data from dilution experiment for dropout analysis.
#' Text file with exported GeneMapper genotypes table.
#' 
#' @docType data
#' @keywords datasets
#' @name set3
#' @method set3 <- import(import.file=system.file("extdata", "set3.txt", package = "strvalidator"))
#' @format ASCII text file
NULL

#' ESX17 example data for dropout analysis.
#' 
#' Reference profiles for source samples.
#' Text file in GeneMapper format.
#' 
#' @docType data
#' @keywords datasets
#' @name ref3
#' @method ref3 <- import(import.file=system.file("extdata", "ref3.txt", package = "strvalidator"))
#' @format ASCII text file
NULL

#' ESX17 example data for dropout analysis.
#' 
#' A slimmed dataset containing data from
#' dilution experiment for dropout analysis (from set3).
#' One sample replicate has lower case sample name (bc9).
#' 
#' @docType data
#' @keywords datasets
#' @name set4
#' @usage data(set4)
#' @format A data frame with 1609 rows and 5 variables
NULL

#' ESX17 example data for dropout analysis.
#' 
#' A slimmed dataset containing reference profiles for source samples in set4.
#' Reference 'A2' has douoble entries for homozygotes.
#' Reference 'F2' has single entries for homozygotes.
#' Reference 'bc' has douoble entries for homozygotes, and lower case sample name. 
#' 
#' @docType data
#' @keywords datasets
#' @name ref4
#' @usage data(ref4)
#' @format A data frame with 98 rows and 3 variables
NULL

#' ESX17 example data for mixture analysis.
#' 
#' A slimmed dataset containing data from
#' mixture experiment for Mx analysis.
#' 
#' @docType data
#' @keywords datasets
#' @name set5
#' @usage data(set5)
#' @format A data frame with 1663 rows and 7 variables
NULL

#' ESX17 example data for mixture analysis.
#' 
#' A slimmed dataset containing the reference profile for the major
#' component  in set5.
#' 
#' @docType data
#' @keywords datasets
#' @name ref51
#' @usage data(ref51)
#' @format A data frame with 34 rows and 3 variables
NULL

#' ESX17 example data for mixture analysis.
#' 
#' A slimmed dataset containing the reference profile for the minor
#' component  in set5.
#' 
#' @docType data
#' @keywords datasets
#' @name ref52
#' @usage data(ref52)
#' @format A data frame with 34 rows and 3 variables
NULL
