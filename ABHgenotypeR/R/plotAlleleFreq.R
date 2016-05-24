#' Plot the  parental allele frequencies along the chromosomes.
#'
#' @param genos The output of readABHgenotypes
#'
#' @return A plot of parental allele frequencies along the chromosomes. If the
#'   output is assigned a name a ggplot2 object is returned for further
#'   manipulation.
#'
#' @examples \dontrun{plotAlleleFreq(genotypes)}
#' \dontrun{p <- plotAlleleFreq(genotypes)}

#' @export
plotAlleleFreq <- function(genos = "genotypes"){

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ####get the number of A,B,H for each marker####
  countAs <- rep(0,length(genos$marker_names))
  for(markerloop in 1:length(genos$marker_names)){

    countAs[markerloop] <- sum(genos$ABHmatrix[,markerloop] == "A", na.rm = TRUE)

  }
  countBs <- rep(0,length(genos$marker_names))
  for(markerloop in 1:length(genos$marker_names)){

    countBs[markerloop] <- sum(genos$ABHmatrix[,markerloop] == "B", na.rm = TRUE)

  }
  countHs <- rep(0,length(genos$marker_names))
  for(markerloop in 1:length(genos$marker_names)){

    countHs[markerloop] <- sum(genos$ABHmatrix[,markerloop] == "H", na.rm = TRUE)

  }
  countNs <- rep(0,length(genos$marker_names))
  for(markerloop in 1:length(genos$marker_names)){

    countNs[markerloop] <- sum(genos$ABHmatrix[,markerloop] == "N", na.rm = TRUE)

  }
  #####

  #####build the dataframe for plotting#####
  genotypeAbs <- cbind("chrom" = genos$chrom,
                       "pos" = genos$pos,
                       "countA" = countAs,
                       "countB" = countBs,
                       "countH" = countHs,
                       "countN" = countNs)
  row.names(genotypeAbs) <- genos$marker_names
  genotypeAbs <- reshape2::melt(as.data.frame(genotypeAbs),id.vars = c("chrom","pos"))

  genotypeRatio <- cbind("chrom" = genos$chrom,
                         "pos" = genos$pos,
                         "ratioA" = 100 / (countAs+countBs+countHs) * countAs,
                         "ratioB" = 100 / (countAs+countBs+countHs) * countBs,
                         "hetero" = 100 / (countAs+countBs+countHs) * countHs)

  row.names(genotypeRatio) <- genos$marker_names
  genotypeRatio <- reshape2::melt(as.data.frame(genotypeRatio),id.vars = c("chrom","pos"),
                        variable.name = "allele_state")
  #####

  #plot the data

  pos <- value <- allele_state <- NULL #appease R cmd check

  ggplot(data = genotypeRatio, aes(x = pos/1000000, y = value, color = allele_state))+
    geom_line(size = 0.75)+
    facet_grid(chrom~.)+
    xlab("physical position (Mb)")+ylab("parental allele frequency (%)")+
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75))+
    scale_color_manual(values = cbbPalette, name = "parental allele",
                      labels = c(genos$nameA,genos$nameB,"hetero"))+

    theme(axis.title = element_text(size = 12, face = "bold"),
          axis.title.y=element_text(vjust = 1.1),
          axis.text = element_text(color = "black", size = 10),
          panel.border = element_rect(fill= NA, colour="grey30"),
          panel.margin = unit(0.2, "lines"),
          axis.ticks = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text.y = element_text(size = 12, angle = 0),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12))
}
