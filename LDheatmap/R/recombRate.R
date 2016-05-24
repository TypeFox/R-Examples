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



#________________Plot recomb Rate track from UCSC genome Browser_______________##
recombRate <- function(minRange, maxRange, chromosome, genome="hg19", vp=viewport(x=0, y=0.99, height=0.04, just=c("left", "top")), view="dense") {
  map_len <- convertX(vp$width, "npc", valueOnly=TRUE)
  Range <- maxRange - minRange
  vp$xscale <- c(minRange, maxRange)
  vp$name <- "recombRate"

  require(rtracklayer)
  require(GenomicRanges)
  session <- browserSession()
  genome(session) <- genome
  query<-ucscTableQuery(session, "recombRate", GRangesForUCSCGenome(genome, chromosome,IRanges(minRange, maxRange)))
  t<-getTable(query)
  if (!dim(t)[1]) {	# we got an empty table
    print ("The genetic region of the data does not correspond to any recombination rate data in the UCSC genome browser")
    return()
  }

  recombRate_title <- textGrob("Recombination Rate from deCODE",
	gp=gpar(fontsize=7, fontfamily="mono"), y=1, just=c("centre", "center"),
 				name="recombRate_title", default.units="native")

  vp_height <- convertX(vp$height, "npc", valueOnly=TRUE)
  rect_height <- 0.5
  if (view == "squish") {
     rect_height <- 0.25
     vp$height <- unit(vp_height*2/3, "npc")
  }
  if (view == "full") {
     rect_height <- 1
     vp$yscale <- c(1+dim(t)[1] ,0)
     recombRate_title <- editGrob(recombRate_title, y=unit(0, "native"))
     vp$height <- unit((1+dim(t)[1])*vp_height/2, "npc")
  }
  if (view == "pack") {
     # Determine the plot line number for each row in t and save it in column t$plot_line
     t$plot_line <- 0
     t$plot_line[1] <- 1
     plot_lines_no <- 1
     if (dim(t)[1] > 1) {	# more than a single gene in table
        text_width <- Range / map_len * convertWidth(grobWidth(
          textGrob(" 8.8 cM/Mb (Avg) ", gp=gpar(fontsize=7, fontfamily="mono"), vp=vp)), "npc", valueOnly=TRUE)
        for (i in 2:dim(t)[1]) {  # for every gene in table starting with the second one
          for (j in 1:plot_lines_no) {
             # compare beggining of region minus text_width to highest end of gene in line j
             if (max(t[1:(i-1),][t[,"plot_line"]==j,"chromEnd"], na.rm = TRUE)
                                         < t[i,"chromStart"] - text_width) {
                t[i, "plot_line"] <- j
                break
             }
          }
          if (!t[i, "plot_line"]) 
             t[i, "plot_line"] <- plot_lines_no <- plot_lines_no + 1
        }
     }
     rect_height <- 1
     vp$yscale <- c(2+plot_lines_no ,0)
     vp$height <- unit((2+plot_lines_no)*vp_height/2, "npc")
  }
     
  recombRate <- gTree(children=gList(recombRate_title), name="recombRate", vp=vp)


  # Read recomb rate for every every row in table t
  for (i in 1:dim(t)[1]) {

       if (view == "dense")       y <- 0  
       else if (view == "squish") y <- rect_height * (i%%2)
       else if (view == "full")   y <- rect_height * i
       else                       y <- rect_height * (1+ t[i, "plot_line"])

       value <- as.integer(t[i,"decodeAvg"]*200)
       color <- as.character(cut(value, breaks=c(0,seq(length=9, from=320, to=600), 10000),
		labels=grey.colors(10)[10:1]))
       rect <- rectGrob(x=max(minRange, t[i, "chromStart"]), default.units="native",
       width=min(maxRange, t[i, "chromEnd"]) - max(minRange, t[i, "chromStart"]),
       y=y, height=rect_height, just=c("left", "bottom"), gp=gpar(col="transparent", fill=color),
	name=paste("rect",i, sep=""))
       recombRate <- addGrob(recombRate , rect)
       if (view == "full" || view == "pack") {
          x <- 1;        if (view == "pack") { x <- i }
          rectText <- textGrob(x=max(minRange, t[x, "chromStart"]), y=y, 
		paste(round(t[i,"decodeAvg"],1)," cM/Mb (Avg) ",sep=""), just=c("right","top"), 
        	gp=gpar(fontsize=7, fontfamily="mono"), default.units="native",
		name=paste("text",i, sep=""))
          recombRate <- addGrob(recombRate , rectText)
       }
  }
  grid.draw(recombRate)
  return(recombRate)
}

