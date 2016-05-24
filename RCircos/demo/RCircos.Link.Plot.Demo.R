#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw chromosome ideogram with padding between 
#       chromosomes, highlights, chromosome names, and link lines. 
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>


    #   Load RCircos library
    #   ============================================
    library(RCircos);


    #   Load human cytoband data and link data
    #   ============================================
    data(RCircos.Link.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   ============================================
    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);

    #   Open the graphic device (here a pdf file)
    #   ============================================
    out.file <- "RCircos.Link.Plot.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();

    #   Draw chromosome ideogram
    #   ============================================
    message("Draw chromosome ideogram ...\n");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Link Plot Demo");

    #   Link lines. Link data has only paired chromosome 
    #   locations in each row and link lines are always 
    #   drawn inside of chromosome ideogram.
    #   ============================================
    message("Add link track ...\n");

    link.data <- RCircos.Link.Data;
    link.colors <- rep("blue", nrow(link.data));
    rows <- seq(1, nrow(link.data), by=5);
    link.colors[rows] <- "red";
    link.data["PlotColor"] <- link.colors;

    track.num <- 2;
    RCircos.Link.Plot(link.data, track.num, FALSE);

    #   Close the graphic device and clear memory
    #   ============================================
    dev.off();
    message("RCircos Link Plot Demo Done!\n");


