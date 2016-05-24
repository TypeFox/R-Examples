#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw chromosome ideogram with padding between  
#       chromosomes, highlights, chromosome names, and heatmap 
#       with gene expression data.
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>


    #   Load RCircos library
    #   =================================================
    library(RCircos);

    #   Load human cytoband data and gene expression data
    #   =================================================
    data(RCircos.Heatmap.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   =================================================
    RCircos.Set.Core.Components(cyto.info, NULL, 5, 0);

    #	Open the graphic device (here a pdf file)
    #   =================================================
    out.file <- "RCircos.Heatmap.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();

    #   Draw chromosome ideogram
    #   =================================================
    message("Draw chromosome ideogram ...");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Heatmap Plot Demo");

    #   Plot five tracks of heatmap
    #   =================================================
    total.track <- 5;
    for(a.track in 1:total.track)
    {
        data.col <- a.track + 4;
        RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, 
                    a.track, "in");
    }

    #   Close the graphic device
    #   =================================================
    dev.off();
    message("RCircos Heatmap demo done!\n");




