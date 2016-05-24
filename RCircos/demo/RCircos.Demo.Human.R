# ______________________________________________________________________
# <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#   This demo draw human chromosome ideogram and data tracks for:
#
#       1.  Connectors
#       2.  Gene lables
#       3.  Heatmap
#       4.  Scatterplot 
#       5.  Line plot
#       6.  Histogram
#       7.  Tile plot
#       8.  Link lines
#       9.  Ribbons
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos package and defined parameters
    #   ===========================================
    library(RCircos);

    #   Load human cytoband data 
    #   ===========================================

    data(UCSC.HG19.Human.CytoBandIdeogram);
    hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #
    #   1.  Chromosome ideogram plot information
    #   2.  x and y coordinates for a circular line and degrees 
    #       of the text rotation at each point
    #   3.  Plot parameters for image plot control
    #  
    #   ======================================================

    RCircos.Set.Core.Components(cyto.info=hg19.cyto, 
        chr.exclude=NULL, tracks.inside=10, tracks.outside=0);

    #   Open the graphic device (here a pdf file)
    #  
    #   ======================================================
    message("Open graphic device and start plot ...");
    pdf(file="RCircos.Demo.Human.pdf", height=8, width=8);

    RCircos.Set.Plot.Area();
    title("RCircos 2D Track Plot with Human Genome");

    #   Draw chromosome ideogram
    #  
    #   ======================================================
    message("Draw chromosome ideogram ...");
    RCircos.Chromosome.Ideogram.Plot();

    #   Connectors and gene names 
    #   ======================================================
    message("Add Gene and connector tracks ...");
    data(RCircos.Gene.Label.Data);

    RCircos.Gene.Connector.Plot(genomic.data=RCircos.Gene.Label.Data, 
                        track.num=1, side="in");
    RCircos.Gene.Name.Plot(gene.data=RCircos.Gene.Label.Data, 
                        name.col=4, track.num=2, side="in");

    #   Heatmap plot.  Since some gene names plotted above are 
    #   wider than one track height, we skip two tracks 
    #   ======================================================
    message("Add heatmap track ...");

    data(RCircos.Heatmap.Data);
    RCircos.Heatmap.Plot(heatmap.data=RCircos.Heatmap.Data, 
            data.col=6, track.num=5, side="in");

    #   Scatterplot
    #   ======================================================
    message("Add scatterplot track ...");

    data(RCircos.Scatter.Data);
    RCircos.Scatter.Plot(scatter.data=RCircos.Scatter.Data, 
            data.col=5, track.num=6, side="in", by.fold=1);

    #   Line plot.
    #   ======================================================
    message("Add line plot track ...");

    data(RCircos.Line.Data);
    RCircos.Line.Plot(line.data=RCircos.Line.Data, data.col=5, 
            track.num=7, side="in");

    #   Histogram plot
    #   ======================================================
    message("Add histogram track ...");

    data(RCircos.Histogram.Data);
    RCircos.Histogram.Plot(hist.data=RCircos.Histogram.Data, 
            data.col=4, track.num=8, side="in");

    #   Tile plot. Note: tile plot data have only chromosome 
    #   locations and each data file is for one track
    #   ======================================================
    message("Add tile track ...");

    data(RCircos.Tile.Data);
    RCircos.Tile.Plot(tile.data=RCircos.Tile.Data, track.num=9, 
                side="in");

    #   Link lines. Link data has only paired chromosome 
    #   locations in each row and link lines are always drawn 
    #   inside of chromosome ideogram.
    #   ======================================================
    message("Add link track ...");

    data(RCircos.Link.Data);
    RCircos.Link.Plot(link.data=RCircos.Link.Data, track.num=11, 
        by.chromosome=FALSE);

    #   Add ribbon link to the center of plot area (link lines).
    #   ======================================================
    message("Add ribbons to link track ...");

    data(RCircos.Ribbon.Data);
    RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, 
            track.num=11, by.chromosome=FALSE, twist=FALSE);

    #   Close the graphic device and clear memory
    #   ======================================================

    dev.off();
    message("R Circos Demo Done ...\n");














