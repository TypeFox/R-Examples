#' EcoSimR
#' EcoSimR is a collection of functions for calculating community metrics and algorithms for randomizing community data for null model analysis. Current modules are included for the analysis of niche overlap, body size overlap, and species co-occurrence. EcoSimR also allows users to define their own functions and algorithms to develop new null models.
#' @name EcoSimR
#' @docType package
NULL

#' MacArthur's (1958) warbler data
#' This data matrix is from MacArthur's classic (1958) paper on the coexistence of 5 species of New England warbler. Each row of the data matrix is a different species of warbler, and each column is one of 16 different subregions of a coniferous tree. Each entry is the percentage of time that each species was observed foraging in a different subregion of the tree (see Figures 2-4 in MacArthur 1958). Zeroes indicate subregions of the tree in which a species was not recorded foraging.
#' @references MacArthur, R.H. 1958. Population ecology of some warblers of northeastern coniferous forests. Ecology 39: 599-699.
#' @name dataMacWarb
#' @docType data
#' @keywords datasets
#' @keywords data
NULL

#' West Indian Finches data
#' This data frame is a binary presence-absence matrix for West Indies finches (Fringillidae). Each row is a different species of finch and each column is one of the 19 major islands in the West Indies. Entries indicate the presence (1) or absence (0) of a species on an island. Data from Gotelli and Abele (1982).
#' @references Gotelli, N.J. and L.G. Abele. 1982. Statistical distributions of West Indian land bird families. Journal of Biogeography 9: 421-435.
#' @name dataWiFinches
#' @docType data
#' @keywords datasets
#' @keywords data
NULL

#' Desert rodent data set
#' This data vector is from Brown's (1975) study of the coexistence of desert rodent species. Each entry is the average adult body mass in grams of six co-occurring species of Sonoran Desert rodents.
#' @references Brown, J.H. 1975. Geographical ecology of desert rodents. p. 314-341 in: Ecology and Evolution of Communities. M.L. Cody and J.M. Diamond (eds.). Harvard University Press, Cambridge.
#' @name dataRodents
#' @docType data
#' @keywords datasets
#' @keywords data
NULL