#' Digitalized Data from a Fluidigm Array
#' 
#' These are the results data from the \code{White} data as measured by the
#' UT digital PCR on Fluidigm 12.765 digital Array. The data were digtilized
#' from a supplementary figure "1471-2164-10-116-S1.pdf" 
#' by White et al. (2009) BMC Genomics
#' 
#' Setup: Experimental details were described be White et al. (2009) BMC 
#' Genomics. The digitalization of the figure was done with imageJ and the 
#' "MicroArray Profile" plugin by Bob Dougherty (rpd@@optinav.com) and 
#' Wayne Rasband.
#'
#' Annotation: See the White et al. (2009) BMC Genomics paper for details.
#' 
#' @name White
#' @docType data
#' @format 
#' \describe{ A dataframe with 9180 rows and 10 columns.
#'  \item{Image_position}{Position of an array in the figure 
#'  1471-2164-10-116-S1.pdf from White et al. (2009) BMC Genomics (e.g., 
#'  11 is the image in the first colum and the first row, 24 is second column 
#'  and fourth image)}
#'  \item{Sample}{is the sample (e.g., "Ace 1:100") as described by White et 
#'                al. (2009) BMC Genomics}
#'  \item{X.1}{Running index for *all* samples}
#'  \item{Index}{Index within an array}
#'  \item{Row}{Row within an array}
#'  \item{Column}{Column within an array}
#'  \item{Area}{is the area that was measured with "MicroArray Profile"}
#'  \item{Min}{is the minimum intensity of an area that was measured with 
#'             "MicroArray Profile"}
#'  \item{Max}{is the maximum intensity of an area that was measured with
#'             "MicroArray Profile"}
#'  \item{Mean}{is the mean intensity of an area that was measured with
#'             "MicroArray Profile"}
#' }
#' @author Stefan Roediger, Michal Burdukiewcz, White et al. (2009) BMC Genomics
#' @references White RA, Blainey PC, Fan HC, Quake SR. Digital PCR provides 
#' sensitive and absolute calibration for high throughput sequencing. 
#' BMC Genomics 2009;10:116. doi:10.1186/1471-2164-10-116.
#'
#' Dougherty B, Rasband W. MicroArray Profile ImageJ Plugin n.d. 
#' http://www.optinav.com/imagej.html (accessed August 20, 2015).
#'
#' @source Data were digitalized from the supplement material (Additional file 
#' 1. dPCR analysis of mock library control.) "1471-2164-10-116-S1.pdf" 
#' by White et al. (2009) BMC Genomics
#' @keywords datasets
#' @examples
#' 
#'str(White)
#'par(mfrow = c(3,3))
#'
# Plot the panels of the arrays (similar to the supplementary figure
# in White et al. (2009) BMC Genomics).
#'
#'White_data <- sapply(unique(White[["Image_position"]]), function(i)
#'		     White[White[["Image_position"]] == i, "Mean"])
#'
#'assays <- sapply(unique(White[["Image_position"]]), function(i)
#'		 unique(White[White[["Image_position"]] == i, "Sample"]))
#'
#'White_adpcr <- create_dpcr(White_data > 115, n = 765, assay = assays, 
#'			   type = "np", adpcr = TRUE)
#'
#'White_k <- colSums(White_data > 115)
#'
#'sapply(2:4, function(i) {
#'    plot_panel(extract_dpcr(White_adpcr, i))
#'
#'    # Create the ECDF of the image scan data to define
#'    # a cut-off for positive and negative partitions
#'    # Plot the ECDF of the image scan data an define a cut-off
#'    plot(ecdf(White_data[, i]), main = paste0("ECDF of Image Scan Data\n", assays[i]),
#'	xlab = "Grey value", ylab = "Density of Grey values")
#'    abline(v = 115, col = 2, cex = 2)
#'    text(80, 0.5, "User defined cut-off", col = 2, cex = 1.5)
#'
#'    # Plot the density of the dPCR experiment
#'    dpcr_density(k = White_k[i], n = 765, bars = TRUE)
#' }
#')
#' 
#'par(mfrow = c(1,1))
#' 
NULL