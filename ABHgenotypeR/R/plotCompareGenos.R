#' Compare to genotype matrices
#'
#' @param genos_1 Output of readABHgenotypes
#' @param genos_2 Output of readABHgenotypes. Note that both genos object need to
#' have identical numbers of marker x individuals.
#' @param markerToPlot A character vector of marker names which appear in the
#'   plot. Defaults to all.
#' @param individualsToPlot A character vector of individual names which appear
#'   in the plot. Defaults to all.
#' @param chromToPlot A character vector of chromosome names which appear in the
#'   plot. Defaults to all.
#' @param CompColors A character vector of length 2 giving the color names or
#'   values to use for differnt and identical markers.
#'   Defaults to black and orange.
#' @param textSize The size of all text elements in the plot. Useful for making a
#'  nice plot. Defaults to 12.
#' @param showMarkerNames Show the marker names along the x axis. This and
#'  showIndividualnames are useful when you display only a few markers and
#'  want them labeled. Defaults to FALSE.
#' @param showIndividualNames Show individual names along the y axis.
#'
#' @return A graphical comparison of genotypes.
#'
#' @examples \dontrun{plotCompareGenos(preImpGenotypes,postImpGenotypes)}
#' \dontrun{#for more examples see plotGenos()}

#' @export
plotCompareGenos <- function(genos_1 = "genotypes_1",
                             genos_2 = "genotypes_2",
                             markerToPlot = "all",
                             individualsToPlot = "all",
                             chromToPlot = "all",
                             CompColors = c("#000000", "#E69F00"),
                             textSize = 12,
                             showMarkerNames = FALSE,
                             showIndividualNames = FALSE) {

  if(showMarkerNames == TRUE) {textX <- element_text(colour = "black", angle = 90)}
  else{textX <- element_blank()}

  if(showIndividualNames == TRUE) {textY <- element_text(colour = "black")}
  else{textY <- element_blank()}

  comp <- genos_1$ABHmatrix == genos_2$ABHmatrix

  comp_df <- as.data.frame(t(comp*1)) #logical -> numerical -> df
  comp_df$chrom <- genos_1$chrom
  comp_df$marker_names <- factor(genos_1$marker_names, levels = unique(genos_1$marker_names))

  comp_df <- reshape2::melt(comp_df,
                            id.vars = c("marker_names","chrom"),
                            variable.name = "individual_names",
                            value.name = "comp")

  if(markerToPlot[1] != "all") comp_df <- comp_df[comp_df$marker_names %in% markerToPlot,]

  if(individualsToPlot[1] != "all") comp_df <- comp_df[comp_df$individual_names %in% individualsToPlot,]

  if(chromToPlot[1] != "all") comp_df <- comp_df[comp_df$chrom %in% chromToPlot,]

  comp_df$individual_names <- factor(comp_df$individual_names,
                                     levels = rev(levels(comp_df$individual_names)))

  marker_names <- individual_names <- comp <- NULL #appease R cmd check

  ggplot(comp_df, aes(x = marker_names, y = individual_names,  fill = factor(comp)))+
    geom_tile()+
    scale_fill_manual(name = "genotype comparison",
                      values = c("0" = CompColors[1], "1" = CompColors[2]),
                      labels = c("different","identical"))+
    facet_grid(.~chrom, scales = "free", space = "free_x")+
    ylab("individuals")+
    xlab("markers")+
    theme(text = element_text(size = textSize),
          axis.text.x = textX,
          axis.text.y = textY,
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom")
}
