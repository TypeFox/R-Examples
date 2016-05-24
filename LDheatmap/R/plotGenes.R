# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

###########################################################################




#_______________________Grab and Plot genes from the UCSC Genome Browser____________##
plotGenes<- function(minRange, maxRange, chromosome, genome="hg19", plot_lines_distance=0.03,
			vp=viewport(x=0, y=0.99, just=c("left", "top")), splice_variants = TRUE, non_coding = TRUE) {
  map_len <- convertX(vp$width, "npc", valueOnly=TRUE)
  Range <- maxRange - minRange
  vp$xscale <- c(minRange, maxRange)
  vp$name <- "transcriptsVP"

  require(rtracklayer)
  require(GenomicRanges)
  session <- browserSession()
  genome(session) <- genome

  query1<-ucscTableQuery(session, "knownGene", GRangesForUCSCGenome(genome, chromosome, IRanges(minRange, maxRange)))
  t<-getTable(query1)
  if (!dim(t)[1]) {	# we got an empty table
    print ("The genetic region of the data does not correspond to any genes in the UCSC genome browser")
    return()
  } 
  t[,"name"] <- as.character(t[,"name"]) # convert from factors to character strings


  # Exclude non-coding genes
  if (!non_coding) {
    ind <- NULL
    for (i in 1:dim(t)[1])
      if (t[i, "cdsStart"] == t[i, "cdsEnd"]
)

        ind <- c(ind, i)
   if (!is.null(ind)) t <- t[-ind,]
  }

  # Exclude splice variants
  if (!splice_variants) {
     # Get canonical splice variant of each gene
     query2 <- ucscTableQuery(session, "knownGene", GRangesForUCSCGenome(genome,chromosome, IRanges(minRange, maxRange)), table="knownCanonical")
     tcanon<-getTable(query2)

     ind<-(t$name %in% as.character(tcanon$transcript))
     t<-subset(t, ind)
  }


  # Get gene names
  t[, "gene_name"] <- ""
  tbl <- "kgXref"
  query2<-ucscTableQuery(session, "knownGene", 
	GRangesForUCSCGenome(genome,chromosome,IRanges(minRange, maxRange)), table=tbl, names=t[,"name"])
  t1<-getTable(query2)
  t1[,"kgID"]       <- as.character(t1[,"kgID"])       # convert factors to character strings
  t1[,"geneSymbol"] <- as.character(t1[,"geneSymbol"]) 
  
  # Put gene names in table t under column "gene_name"
  for (i in 1:dim(t)[1]) {
    gene_name <- t1[t1[,"kgID"]== t[i,"name"],"geneSymbol"]
    if (length(gene_name) != 0) t[i,"gene_name"] <- gene_name
  }

  # Get gene colors
  if ("kgColor" %in% tableNames(query1)) {
    query3<-ucscTableQuery(session, "knownGene", 
	GRangesForUCSCGenome(genome,chromosome,IRanges(minRange, maxRange)), table="kgColor", 
	names=t[,"name"])
    color_tbl<-getTable(query3)
  }

  # Determine the plot line number for each row in t and save it in column t$plot_line
  t$plot_line <- 0
  t$plot_line[1] <- 1
  plot_lines_no <- 1
  if (dim(t)[1] > 1) {	# more than a single gene in table
  for (i in 2:dim(t)[1]) {  # for every gene in table starting with the second one
    gene_name_width <- Range / map_len * convertWidth(grobWidth(
      textGrob(paste(t[1,"gene_name"], "   "), gp=gpar(fontsize=7))), "npc", valueOnly=TRUE)
    for (j in 1:plot_lines_no) {
       # compare beggining of gene minus gene_name_width to highest end of gene in line j
       if (max(t[1:(i-1),][t[,"plot_line"]==j,"txEnd"], na.rm = TRUE)
                                         < t[i,"txStart"] - gene_name_width) {
          t[i, "plot_line"] <- j
          break
       }
    }
    if (!t[i, "plot_line"]) 
       t[i, "plot_line"] <- plot_lines_no <- plot_lines_no + 1
  }}


  height=plot_lines_no * plot_lines_distance
  vp$height <- unit(height, "npc")


  # Create gene plot title
  gene_plot_title <- textGrob("UCSC Genes Based on RefSeq, UniProt, GenBank, CCDS and Comparative Genomics",
	gp=gpar(fontsize=7, fontfamily="mono"), y=1, just=c("centre", "center"),
 		name="gene_plot_title", default.units="native")

   Transcripts <- gTree(children=gList(gene_plot_title), name="transcripts", vp=vp)


  # Plot every row in table t
  for (i in 1:dim(t)[1]) {
     tx_upArrows <- tx_downArrows <- tx_exons <- tx_cds <- tx_leftArrow <- tx_rightArrow <- NULL
     y = 1 - t[i, "plot_line"] / plot_lines_no
     txStart <- max(minRange, t[i,"txStart"])
     txEnd   <- min(maxRange, t[i, "txEnd"])
     cdsStart <- max(minRange, t[i,"cdsStart"])
     cdsEnd <- min(maxRange, t[i,"cdsEnd"])
     exonStarts <- as.numeric(strsplit(as.character(t[i,"exonStarts"]), ",")[[1]])
     exonStarts[exonStarts < minRange] <- minRange
     exonStarts[exonStarts > maxRange] <- maxRange
     exonEnds <- as.numeric(strsplit(as.character(t[i,"exonEnds"]), ",")[[1]])
     exonEnds[exonEnds < minRange] <- minRange
     exonEnds[exonEnds > maxRange] <- maxRange

     # Print gene name
     tx_name <- textGrob(x=txStart, y=y, paste(t[i,"gene_name"]," ",sep=""), just=c("right","center"), 
        	gp=gpar(fontsize=7, fontfamily="mono"), default.units="native", name="gene_name")

     # Plot UCSC transcription region
     tx_region <- linesGrob(x=c(txStart, txEnd), y=c(y,y), gp=gpar(lwd=1, lineend="butt"),
 		default.units="native", name="tx_region")

     # Plot arrow heads along the transcription region
     distance_between_arrows <- Range * (10/6) /
		 convertX(unit(map_len,"npc"), "millimeters", valueOnly=TRUE)
     no_arrows <- floor(1 / distance_between_arrows * (txEnd - txStart)) - 1

     if (no_arrows > 0) {
     arrow_heads_x0 <-unit(txStart + 1:no_arrows*distance_between_arrows, "native")
     arrow_heads_x1 <-unit(txStart + 1:no_arrows*distance_between_arrows, "native") + unit(0.5, "millimeters")
     if ( t[i,"strand"] == "+") { 
        change_y0 = unit(0.5, "millimeters"); change_y1 <- unit(0, "native") } 
     else { 
        change_y0 <- unit(0, "native"); change_y1 <- unit(0.5, "millimeters") }
     tx_upArrows <- segmentsGrob( x0=arrow_heads_x0, x1=arrow_heads_x1,
        y0=unit(y, "native")+change_y0, y1=unit(y, "native")+change_y1, 
        name="tx_upArrows")
     tx_downArrows <- segmentsGrob( x0=arrow_heads_x0, x1=arrow_heads_x1,
        y0=unit(y, "native")-change_y0, y1=unit(y, "native")-change_y1, 
        name="tx_downArrows")
     }

     # Plot UCSC exons
     tx_exons <- segmentsGrob(x0=exonStarts, x1=exonEnds, y0=y, y1=y, gp=gpar(lwd=5,lineend="butt"), 					default.units = "native", name="tx_exons")

     # Calculate where to draw coding region
     cds0 <- cds1 <- NULL
     for (j in 1:length(exonStarts)) {
        if (cdsStart < exonEnds[j] & cdsEnd > exonStarts[j]) {
           cds0 <- c(cds0, max(exonStarts[j], cdsStart))
           cds1 <- c(cds1, min(exonEnds[j], cdsEnd))
     }}

     # Plot UCSC coding region
     if(!is.null(cds0))		# draw only if coding region is within plot region
        tx_cds <- segmentsGrob(x0=cds0, x1=cds1, y0=y, y1=y, 
		gp=gpar(lwd=10,lineend="butt"), default.units="native", name="tx_cds")

     # Plot two white arrow heads pointing to the left if necessary
     if (t[i,"txStart"] < minRange)
         tx_leftArrow <- polygonGrob(x=	unit.c(unit(minRange, "native")+unit(1, "millimeters"), 
			unit(minRange, "native")+unit(2, "millimeters"), 
			unit(minRange, "native")+unit(2, "millimeters"),
			unit(minRange, "native")+unit(2, "millimeters"),
			unit(minRange, "native")+unit(3, "millimeters"),
			unit(minRange, "native")+unit(3, "millimeters")), 
                      y=	rep(unit.c(unit(y, "native"), 
			unit(y, "native")+unit(1, "millimeters"), 
			unit(y, "native")-unit(1, "millimeters")),2), 
                      rep(1:2, each=3), gp=gpar(fill="white"), name="tx_leftArrow")

     # Plot two white arrow heads pointing to the right if necessary
     if (t[i,"txEnd"] > maxRange)
         tx_rightArrow <- polygonGrob(x= unit.c(unit(maxRange, "native")-unit(1, "millimeters"), 
			unit(maxRange, "native")-unit(2, "millimeters"), 
			unit(maxRange, "native")-unit(2, "millimeters"),
			unit(maxRange, "native")-unit(2, "millimeters"),
			unit(maxRange, "native")-unit(3, "millimeters"),
			unit(maxRange, "native")-unit(3, "millimeters")), 
                      y=	rep(unit.c(unit(y, "native"), 
			unit(y, "native")+unit(1, "millimeters"), 
			unit(y, "native")-unit(1, "millimeters")),2), 
                      rep(1:2, each=3), gp=gpar(fill="white"), name="tx_rightArrow")

     # Determine gene color
     gene_color <- "grey50"
     if (exists("color_tbl")) {
        r <- color_tbl[color_tbl[,"kgID"]== t[i,"name"],"r"]
        g <- color_tbl[color_tbl[,"kgID"]== t[i,"name"],"g"]
        b <- color_tbl[color_tbl[,"kgID"]== t[i,"name"],"b"]
        gene_color <- rgb(r,g,b, maxColorValue = 225)
     }

     # Wrap everything in a grob and draw it
     transcript <- gTree(children =
        gList(tx_name, tx_region, tx_upArrows, tx_downArrows, tx_exons, tx_cds,tx_leftArrow, 
	 tx_rightArrow), gp=gpar(col=gene_color), name=t[i, "name"])
     Transcripts <- addGrob(Transcripts, transcript)
 }

  grid.draw(Transcripts)
  return(Transcripts)
}

