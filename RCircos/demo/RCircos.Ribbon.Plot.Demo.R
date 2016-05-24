#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw ribbon links between chromosomes in the center
#       of chromosome ideogram
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>


    #   Load RCircos library
    #   ===========================================
    library(RCircos);

    #   Load human cytoband data and link data
    #   ===========================================
    data(RCircos.Ribbon.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   ===========================================
    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


    #   Open the graphic device (here a pdf file)
    #   ===========================================
    out.file <- "RCircos.Ribbon.Plot.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();

    #   Draw chromosome ideogram
    #   ===========================================
    message("Draw chromosome ideogram ...\n");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Ribbon Plot Demo");


    #   Link lines. Link data has only paired chromosome 
    #   locations in each row and link lines are always 
    #   drawn inside of chromosome ideogram.
    #   ================================================
    message("Add link track ...");

    ribbon.data <- RCircos.Ribbon.Data;
    ribbon.data["PlotColor"] <- c("red", "blue", "green", "cyan");

    track.num <- 11;
    RCircos.Ribbon.Plot(ribbon.data, track.num, FALSE);


    #   Close the graphic device 
    #   ===========================================
    dev.off();
    message("RCircos Ribbon Plot Demo Done!\n");





