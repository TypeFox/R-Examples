#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw chromosome ideogram with padding between 
#       chromosomes, highlights, chromosome names, and scatterplot. 
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>


    #   Load RCircos library
    #   ===============================================
    library(RCircos);

    #   Load human cytoband data and scatterplot data
    #   ===============================================
    data(RCircos.Scatter.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   ===============================================
    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);

    #   Open the graphic device (here a pdf file)
    #   ===============================================
    out.file <- "RCircos.Scatter.Plot.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();

    #	Draw chromosome ideogram
    #   ===============================================
   message("Draw chromosome ideogram ...");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Scatter Plot Demo");

    #   Scatterplot 
    #   ===============================================
    scatter.data <- RCircos.Scatter.Data;
    scatter.colors <- rep("cyan", nrow(scatter.data));
    scatter.colors[which(scatter.data$seg.mean>=2)] <- "red";
    scatter.colors[which(scatter.data$seg.mean<=-2)] <- "blue";
    scatter.data["PlotColor"] <- scatter.colors;

    data.col <- 5;
    track.num <- 6; 
    RCircos.Scatter.Plot(scatter.data, data.col, track.num, "in", 1);

    #   Close the graphic device
    #   ===============================================
    dev.off();
    message("RCircos Scatter Plot Demo Done!");

