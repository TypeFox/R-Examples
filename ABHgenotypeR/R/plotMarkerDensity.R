#' Plot the marker density along the chromosomes.
#'
#' @param genos The output of readABHgenotypes
#'
#' @return A plot of marker densities along the chromosomes. If the output is
#'   assigned a name a ggplot2 object is returned for further manipulation.
#'
#' @examples \dontrun{plotMarkerDensity(genotypes)}
#' \dontrun{p <- plotMarkerDensity(genotypes)}

#' @export
plotMarkerDensity <- function(genos = "genotypes"){

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  SNPdistr <- cbind.data.frame("chrom" = genos$chrom,
                               "pos" = genos$pos)

  pos <- NULL #appease R cmd check

  ggplot(data = SNPdistr, aes(x = pos/1000000))+
    stat_bin(binwidth = 1, drop = TRUE, geom = "line", size = 0.75)+
    labs(x = expression(bold(physical~position~(Mb))),
         y = expression(bold(site~density~(no.~of~sites~Mb^{-1}))))+
    facet_grid(chrom~.)+
    scale_colour_manual(values=cbbPalette, name = "max missing")+
    theme(axis.title = element_text(size = 12), # bold is wrapped in expression()
          axis.text = element_text(color = "black", size = 10),
          panel.border = element_rect(fill= NA, colour="grey30"),
          axis.ticks = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text.y = element_text(size = 12, angle = 0),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12))
}
