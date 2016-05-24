#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw chromosome ideogram with padding area between 
#       chromosomes, highlights, and chromosome names, label gene 
#       names in the second data track and put connectors between 
#       chromosome ideogram and gene names.
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos package
    #   ============================================
    library(RCircos);

    #   Load gene label data and human cytoband data 
    #   ============================================
    data(RCircos.Gene.Label.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   ============================================
    RCircos.Set.Core.Components(cyto.info, NULL, 5, 5);

    #   Ready to make Circos plot image. Open graphic 
    #   device (pdf file)
    #   ============================================
    out.file <- "RCircos.Gene.Connector.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();
    title("RCircos Gene and Connector Plot Demo");

    #   Draw chromosome ideogram
    #   ============================================
    message("Draw chromosome ideogram ...");
    RCircos.Chromosome.Ideogram.Plot();

    #   Plot connectors and gene names inside of 
    #   chromosome ideogram. Connectors are in first 
    #   track and gene names in the second track. 
    #   ============================================
    message("Add Gene and connector tracks ...");
    data(RCircos.Gene.Label.Data);

    direction <- "in";
    track.num <- 1;
    RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, 
            track.num, direction);

    gene.data <- RCircos.Gene.Label.Data;
    gene.colors <- rep("black", nrow(gene.data))
    gene.colors[which(gene.data$Gene=="TP53")] <- "red";
    gene.colors[which(gene.data$Gene=="BRCA2")] <- "red";
    gene.colors[which(gene.data$Gene=="RB1")] <- "red";
    gene.colors[which(gene.data$Gene=="JAK1")] <- "blue";
    gene.colors[which(gene.data$Gene=="JAK2")] <- "blue";

    gene.data["PlotColor"] <- gene.colors;

    name.col <- 4;
    track.num <- 2;
    RCircos.Gene.Name.Plot(gene.data, name.col, track.num, direction);

    #   Plot connectors and gene names outside of 
    #   chromosome ideogram. Connectors are in first 
    #   track and gene names in the second track. 
    #   ============================================
    direction <- "out";
    track.num <- 1;
    RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, 
            track.num, direction);
    track.num <- 2;
    RCircos.Gene.Name.Plot(gene.data, name.col, track.num, 
            direction);

    #   Close the graphic device and clear memory
    #   ============================================
    dev.off();
    message("R Circos Demo Done ...\n");


