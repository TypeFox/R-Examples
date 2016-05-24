#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw chromosome ideogram with padding between 
#       chromosomes, highlights, chromosome names, and histogram. 
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos package and defined parameters
    #   ==================================================
    library(RCircos);

    #   Load human cytoband data and gene expression data
    #   ==================================================
    data(RCircos.Histogram.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   ================================================
    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);

    #   Open the graphic device (here a pdf file)
    #   ==================================================
    out.file <- "RCircos.Histogram.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();

    #   Draw chromosome ideogram
    #   ==================================================
    message("Draw chromosome ideogram ...");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Histogram Plot Demo");

    #   Plot histogram Inside of chromosome ideogram
    #   ==================================================
    hist.data <- RCircos.Histogram.Data;
    hist.colors <- rep("blue", nrow(hist.data))
    rows <- which(hist.data$Data>0.4)
    hist.colors[rows] <- "red";
    hist.data["PlotColor"] <- hist.colors;

    data.col <- 4; 
    track.num <- 1;
    RCircos.Histogram.Plot(hist.data, data.col, track.num, "in");

    #   Close the graphic device and clear memory
    #   ==================================================
    dev.off();
    print("RCircos Histogram Demo Done!\n");

