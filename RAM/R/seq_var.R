seq_var <- function (taxon=NULL, region="ITS", align=NULL, file.align=NULL, 
                     file.align.format="fasta", outgroup.name=NULL, 
                     taxa.sep=":", col.total=5, col.genus=1, col.genus.syn=3,
                     col.species=4, col.strain=5,
                     box.cex.axis.text.x=6, box.cex.axis.text.y=6, 
                     box.cex.xlab=8, box.cex.ylab=8, 
                     intra.fill="purple", inter.fill="orange",
                     den.cex.axis.text.x=6, den.cex.axis.text.y=6,
                     den.cex.xlab=6, den.cex.ylab=6, frame.col="blue",
                     table.cex=8, main=FALSE, main.cex=14,
                     file=NULL, width=8, height=8) {
  #library(ape)
  #library(phangorn)
  #library(gdata)
  
  if (!requireNamespace("gridExtra") ) {
  #if (!require("gridExtra") ) {
    stop("R package gridExtra is required by this function")
  } 
 
  sys <- switch(Sys.info()[['sysname']],
  Windows= {print("I'm a Windows PC.")},
  Linux  = {print("I'm a penguin (Linux).")},
  Darwin = {print("I'm a Mac.")})
  if ( sys == "I'm a Windows PC." ) {
    stop("Sorry, I am not functional. This R function works on Linux and Mac operating systems")
  } else {
    print("OK, this is your lucky day, I am working")
  }

   save <- !is.null(file)

  # taxon
  if ( is.null(taxon) ) {
    stop("please provide the taxonomic information of the alignment")
  }

  # Provide align or Read fasta file into "align"
  if ( is.null(align) & is.null(file.align) ) {
    stop("must provide alignment variable or a file of alignment in fasta format")
  }

  if ( !is.null(align) & !is.null(file.align) ) {
    stop ("Can only provide one of the following: alignment variable or a file of alignment in fasta format")
  }

  if ( is.null(align) ) {
    #align_all <- read.dna(FILE_NAME, format="fasta")
    align_all <- ape::read.dna(file.align, format=file.align.format)
  } else {
    align_all <- align
  }
  #return(align_all)
  # exclude outgroup
  
  if (!is.null(outgroup.name) ) {
    if (length(outgroup.name)==1L) {
       rm <- paste0(":", outgroup.name, ":")
    } else {
       rm <- paste(":", outgroup.name, ":", collapse="|", sep="")
    }
    align1 <- align_all[!grepl(rm,rownames(align_all)), ]
  } else {
    align1 <- align_all
  }
  #return(align1)
  # raw distance is p-distance with substitution -> d: transition + transversion
  system.time(dm <- ape::dist.dna(align1, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE) )

  dm[upper.tri(dm, diag = TRUE)] <- NA
  if ( !requireNamespace("reshape2") ) {
    stop("R package reshape2 is required for this function")
  }
  dm.m <- reshape2::melt(dm)
  # remove na
  dm.m <- dm.m[!is.na(dm.m$value),]
  #dim(dm.m)
  # check whether the names are factor
  #str(dm.m)
  # convert factor to character before split the columns
  for (i in 1:ncol(dm.m) ) {
    if (is.factor(dm.m[,i])) {
       dm.m[,i] <- as.character(dm.m[,i])
    }
  }
  #str(dm.m)
  #return(dm.m)
  # write dm.m out
  write.csv(dm.m, file="temp.dm.m.csv", quote=FALSE)
  
  # split the new columns
  #system.time(dm.m.split <- col.splitup(dm.m.split, col="Var1_split", sep=":", drop=TRUE))
  #system.time(dm.m.split2 <- col.splitup(dm.m.split, col="Var2_split", sep=":", drop=TRUE))
  system("tr < temp.dm.m.csv -s ':' | tr ':' ',' > temp.dm.m.split.csv")
    #return(dm.m)
  dm.m.split <- suppressWarnings(RAM::fread.meta("temp.dm.m.split.csv"))


  # renames columns
  num_col <- 2*col.total+1
  if ( ncol(dm.m.split) != num_col ) {
    stop ("Please make sure the sequence ID of your alignment is something similar to the following format: genus_name:seq_identifier:genus_synonym:species:strain")
  } else { 
    names(dm.m.split)[c(col.genus, col.genus.syn, col.species, col.strain, col.genus + col.total, col.genus.syn + col.total, col.species + col.total, col.strain + col.total, ncol(dm.m.split))] <- c("seq1_genus", "seq1_genus_syn", "seq1_species", "seq1_strain", "seq2_genus", "seq2_genus_syn", "seq2_species", "seq2_strain", "distance")
  # head(dm.m.split)
  write.csv(dm.m.split, file=paste0(taxon, ".dm.m.split.csv"), quote=FALSE)
  }
  
  #return(dm.m.split)
  # Makes Intraspecific subset for Full dataset
  #############################################################################
  # note, by changing symbol to not equal to (!=) you get interspecific data 
  #############################################################################
  # define columns for variation calculation
    # exclude sp_ (undescribed species)
  dim(dm.m.split)
  dm.m.split <- dm.m.split[ dm.m.split[, col.species] != "sp_" &   dm.m.split[, col.species + col.total] != "sp_", ]
  dim(dm.m.split)
  
  # check whether the distance column is numeric 
  if ( !is.numeric(dm.m.split$distance) ) {
    stop ( "the distance column is not numeric, it might because the sequence ID of the alignment is not in the right format, which should be something similar to : genus_name:seq_identifier:genus_synonym:species:strain")
  }

  # intraspecific 
  intra <- dm.m.split[ which((dm.m.split[,"seq1_genus_syn"] == dm.m.split[,"seq2_genus_syn"]) & (dm.m.split[,"seq1_species"] == dm.m.split[,"seq2_species"]) ),]
  dim(intra)
  # Order the data by high level classification and in reverse order of distance
  intra_ord <- intra[order(-as.numeric(intra[,"distance"]), intra[,"seq1_genus_syn"], intra[,"seq1_species"], intra[,"seq1_strain"]), ]
  intra_ord[,"seq1_species"] <- paste( intra_ord[,"seq1_genus_syn"],  intra_ord[,"seq1_species"], sep=" ")
  intra_ord[,"seq2_species"] <- paste( intra_ord[,"seq2_genus_syn"],  intra_ord[,"seq2_species"], sep=" ")

  # interspecific
  inter <- dm.m.split[ which((dm.m.split[,"seq1_genus_syn"] == dm.m.split[,"seq2_genus_syn"])  & (dm.m.split[,"seq1_species"] != dm.m.split[,"seq2_species"]) ),]
  dim(inter)
  # Order the data by high level classification and in reverse order of distance
  inter_ord <- inter[order(-as.numeric(inter[,"distance"]), inter[,"seq1_genus_syn"], inter[,"seq1_species"], inter[,"seq1_strain"]), ]
  inter_ord[,"seq1_species"] <- paste( inter_ord[,"seq1_genus_syn"],  inter_ord[,"seq1_species"], sep=" ")
  inter_ord[,"seq2_species"] <- paste( inter_ord[,"seq2_genus_syn"],  inter_ord[,"seq2_species"], sep=" ")
  
  #return(list(intra=intra, inter=inter))
  # convert summary(intra$distance) to dataframe
  intra.sum<-as.data.frame(t(round(summary(intra_ord$distance), 4)))
  intra.sum <- intra.sum[, -1]
  rownames(intra.sum) <- intra.sum$Var2
  intra.sum <- intra.sum[, -1, drop=FALSE]
  names(intra.sum) <- "Intraspecific"

  # convert summary(inter$distance) to dataframe
  inter.sum<-as.data.frame(t(round(summary(inter_ord$distance), 4)))
  inter.sum <- inter.sum[, -1]
  rownames(inter.sum) <- inter.sum$Var2
  inter.sum <- inter.sum[, -1, drop=FALSE]
  names(inter.sum) <- "Interspecific"

    # library(scales)
  if (!requireNamespace("scales") ) {
     stop("R package scales is required for this function")
  }
  
  # reorder by distance
  #suppressWarnings(intra_ord$seq1_species <- factor(intra_ord$seq1_species , levels = intra_ord[order(intra_ord$distance, decreasing=TRUE), "seq1_species"]))
  #suppressWarnings(intra_ord$seq1_species <- factor(intra_ord$seq1_species , levels = rev(levels(factor(intra_ord$seq1_species)))))
  #suppressWarnings(inter_ord$seq1_species <- factor(inter_ord$seq1_species , levels = inter_ord[order(inter_ord$distance, decreasing=TRUE), "seq1_species"]))
  #suppressWarnings(inter_ord$seq1_species <- factor(inter_ord$seq1_species , levels = rev(levels(factor(inter_ord$seq1_species)))))
  # reorder by median
  ord_intra <- order(as.numeric(by(intra_ord$distance, intra_ord$seq1_species, median)))
  intra_ord$seq1_species <- ordered(intra_ord$seq1_species, levels=levels(factor(intra_ord$seq1_species))[ord_intra])
  
  ord_inter <- order(as.numeric(by(inter_ord$distance, inter_ord$seq1_species, median)))
  inter_ord$seq1_species <- ordered(inter_ord$seq1_species, levels=levels(factor(inter_ord$seq1_species))[ord_inter])  
  
  #return(list(intra=intra_ord, inter=inter_ord))
  # Plot intra-specific variation
  # plot boxplot

  # To handle R CMD check 'no visible binding for global variable' notes, 
  #1) NULLing. e.g. x<-NULL; 
  #2) use with block, e.g. with(data, ggplot(...))
  #3) use aes_string
  ..count..<-NULL
  distance <- NULL
  seq1_species <- NULL

  bp_intra2 <- ggplot(data=intra_ord, aes(y = distance, x = seq1_species)) + 
      geom_boxplot(outlier.size = 1, fill=intra.fill) + coord_flip() +
      xlab("species") + ylab("pairwise distance") + 
      #ggtitle(paste0("Intraspecific variation of ", taxon, " ", region, " region")) + 
      theme(axis.text.y = element_text(size = box.cex.axis.text.y, hjust = 1, face = 'italic'), axis.text.x = element_text(size = box.cex.axis.text.x, face="italic"), axis.title.x=element_text(size=box.cex.xlab,face="bold"), axis.title.y=element_text(size=box.cex.ylab,face="bold") ) 
      #scale_y_reverse()

  # plot density
  bp_intra_density <- ggplot(intra_ord, aes(x=distance)) + 
      geom_density(aes(y=..count../sum(..count..)), alpha=0.2, , fill=intra.fill) +
      coord_flip() + scale_y_continuous(labels = scales::percent_format()) +
      xlab("pairwise distance") + ylab("percent count") +
      theme(axis.text.y = element_text(size = den.cex.axis.text.y, hjust = 1, face = 'italic'), axis.text.x = element_text(size = den.cex.axis.text.x), axis.title.x=element_text(size=den.cex.xlab,face="bold"), axis.title.y=element_text(size=den.cex.ylab,face="bold")) 

  # Plot inter-specific variation
  # plot boxplot
  bp_inter2 <- ggplot(data=inter_ord, aes(y = distance, x = seq1_species)) + 
      geom_boxplot(outlier.size = 1, fill=inter.fill) + coord_flip() +
      xlab("species") + ylab("pairwise distance") + 
      theme(axis.text.y = element_text(size = box.cex.axis.text.y, hjust = 1, face = 'italic'), axis.text.x = element_text(size = box.cex.axis.text.x, face="italic"), axis.title.x=element_text(size=box.cex.xlab,face="bold"), axis.title.y=element_text(size=box.cex.ylab,face="bold")) 
      #scale_y_reverse()

  # plot density
  bp_inter_density <- ggplot(inter_ord, aes(x=distance)) + 
      geom_density(aes(y=..count../sum(..count..)), alpha=0.2, , fill=inter.fill) +
      coord_flip() + scale_y_continuous(labels = scales::percent_format()) +
      xlab("pairwise distance") + 
      ylab("percent count") +
      theme(axis.text.y = element_text(size = den.cex.axis.text.y, hjust = 1, face = 'italic'), axis.text.x = element_text(size = den.cex.axis.text.x), axis.title.x=element_text(size=den.cex.xlab,face="bold"), axis.title.y=element_text(size=den.cex.ylab,face="bold"))  


  if ( main )  {
     title=textGrob(paste("Intra (left) and inter (right) -specific variation of ", taxon, " ", region, " region", sep=""), gp=gpar(fontsize=main.cex, fontface="bold"), just="top", vjust=1)
  } else if ( is.character(main) )  {
     title=textGrob(main, gp=gpar(fontsize=main.cex, fontface="bold"), just="top")
  } else {
     title=textGrob("")
  }

  # table theme
  t_theme <- getFromNamespace("ttheme_default", "gridExtra")(core=list(fg_params=list(fontsize=table.cex)), colhead=list(fg_params=list(fontsize=table.cex, fontface="bold")))

  p1 <- arrangeGrob(bp_intra2)
  intra.sum<-cbind(rownames(intra.sum), intra.sum)
  rownames(intra.sum) <-NULL
  names(intra.sum)[1]<-""
  p2 <- arrangeGrob(tableGrob(intra.sum, rows=NULL, cols=colnames(intra.sum), t_theme), bp_intra_density, ncol=1)
  p3 <- arrangeGrob(grobTree(arrangeGrob(p1, p2, ncol=2, widths=c(5, 2), clip=TRUE), rectGrob(.5, .5, width=unit(.99, "npc"), height=unit(0.99, "npc"), gp=gpar(lwd=2, fill=NA, col=frame.col, alpha=0.2))))

  p4 <- arrangeGrob(bp_inter2)
  inter.sum<-cbind(rownames(inter.sum), inter.sum)
  rownames(inter.sum) <-NULL
  names(inter.sum)[1]<-""
  p5 <- arrangeGrob(tableGrob(inter.sum, rows=NULL, cols=colnames(inter.sum), t_theme), bp_inter_density, ncol=1)
  p6 <- arrangeGrob(grobTree(arrangeGrob(p4, p5, ncol=2, widths=c(5, 2), clip=TRUE), rectGrob(.5, .5, width=unit(.99, "npc"), height=unit(0.99, "npc"), gp=gpar(lwd=2, fill=NA, col=frame.col, alpha=0.2))))

  if ( main | is.character(main) ) {
    p7 <- arrangeGrob(title, arrangeGrob(p3, p6, ncol=2, widths=c(1,1)), nrow=2, heights=c(1, 10), clip=TRUE)
  } else {
    p7 <- arrangeGrob(p3, p6, ncol=2, widths=c(1,1), clip=TRUE)
  }

  if (save) {
    grid.newpage();
    pdf(file=file, width=width, height=height)
    grid.draw(p7)
    dev.off()
    grid.newpage();
    grid.draw(p7);
  } else {
    grid.newpage();
    grid.draw(p7)
  }
  system("rm temp.*.csv")

}


