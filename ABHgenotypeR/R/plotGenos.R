#' Plot graphical genotypes.
#'
#' @param genos The output of readABHgenotypes
#' @param markerToPlot A character vector of marker names which appear in the
#'   plot. Defaults to all.
#' @param individualsToPlot A character vector of individual names which appear
#'   in the plot. Defaults to all.
#' @param chromToPlot A character vector of chromosome names which appear in the
#'   plot. Defaults to all.
#' @param alleleColors A character vector of length 4 giving the color names or
#'   values to use for the A,B,H and n.d genotypes. Defaults to orange, blue,
#'   green and black.
#' @param textSize The size of all text elements in the plot. Useful for making a
#'  nice plot. Defaults to 12.
#' @param showMarkerNames Show the marker names along the x axis. This and
#'  showIndividualNames are useful when you display only a few markers and
#'  want them labeled. Defaults to FALSE.
#' @param showIndividualNames Show individual names along the y axis.
#'
#' @return Graphical genotypes.
#'
#' @examples \dontrun{plotGenos(genotypes)}
#' markerNames <- c("marker1", "marker2", "marker3")
#' individualNames <- c("F2_100", "F2_101", "F2_102", "F2_103")
#' someColors <- c("black", "red", "gold", "white")
#' \dontrun{plotgenos(genotypes, markerNames, individualNames, 1:3, someColors)}
#'
#' \dontrun{p <- plotGenos(genotypes)}

#' @export
plotGenos <- function(genos = "genotypes",
                      markerToPlot = "all",
                      individualsToPlot = "all",
                      chromToPlot = "all",
                      alleleColors = c("#56B4E9","#E69F00",
                                       "#009E73", "#000000"),
                      textSize = 12,
                      showMarkerNames = FALSE,
                      showIndividualNames = FALSE) {

  if(showMarkerNames == TRUE) {textX <- element_text(colour = "black", angle = 90)}
  else{textX <- element_blank()}

  if(showIndividualNames == TRUE) {textY <- element_text(colour = "black")}
  else{textY <- element_blank()}

  ggt <- data.frame(t(genos$ABHmatrix), stringsAsFactors = FALSE, check.names = FALSE)

  ggt$chrom <- genos$chrom
  ggt$index <- 1:length(genos$chrom)
  ggt$marker_names <- factor(genos$marker_names, levels = unique(genos$marker_names))

  ggt <- reshape2::melt(ggt,
                        id.vars = c("chrom", "index","marker_names"),
                        variable.name = "individual_names", value.name = "allele")

  if(markerToPlot[1] != "all") ggt <- ggt[ggt$marker_names %in% markerToPlot,]

  if(individualsToPlot[1] != "all") ggt <- ggt[ggt$individual_names %in% individualsToPlot,]

  if(chromToPlot[1] != "all") ggt <- ggt[ggt$chrom %in% chromToPlot,]

  marker_names <- individual_names <- allele <- NULL #appease R cmd check

  ggt$individual_names <- factor(ggt$individual_names,
                                     levels = rev(levels(ggt$individual_names)))

  ggplot(ggt, aes(x = marker_names, y = individual_names,
                  fill = allele))+
    geom_tile()+
    scale_fill_manual(name = "genotypes",
                      values = c("A" = alleleColors[1], "B" = alleleColors[2],
                                 "H" = alleleColors[3], "N" = alleleColors[4]),
                      labels = c(genos$nameA, genos$nameB, "hetero", "n.d."))+
    facet_grid(.~chrom, scales = "free", space = "free_x")+
    xlab("marker")+
    ylab("individuals")+
    theme(text = element_text(size = textSize),
          axis.text.x = textX,
          axis.text.y = textY,
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position = "bottom",
          panel.margin = unit(0, "lines"))

}
