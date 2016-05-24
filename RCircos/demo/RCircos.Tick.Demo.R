#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
#
#       This demo draw human chromosome ideogram with ticks. Since in 
#       most cases ticks are not necessary for data presentation, so
#       we add ticks to the existing ideogram instead of plotting an 
#       ideogram with ticks.
#
#       Set 20 tracks inside of chromosome ideogram to get enough space 
#       for all chromosomes with ticks in 5MB interval
#   ______________________________________________________________________
#   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>


    #   Load RCircos package
    #   ===========================================
    library(RCircos);

    #   Load human cytoband data 
    #   ===========================================
    data(UCSC.HG19.Human.CytoBandIdeogram);
    hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;

    #   Setup RCircos core components:
    #   ===========================================
    RCircos.Set.Core.Components(cyto.info=hg19.cyto, 
            chr.exclude=NULL, tracks.inside=20, 
            tracks.outside=5);

    #   Open the graphic device (here a pdfimage file)
    #   ==============================================
    pdf("RCircos.Tick.Demo.pdf", height=8, width=8);
    RCircos.Set.Plot.Area();

    #   Draw chromosome ideogram then add ticks with 
    #   5 mb interval
    #   ===========================================
    RCircos.Chromosome.Ideogram.Plot();
    RCircos.Add.Ideogram.Tick(5);
    title("RCircos Tick Demo")

    #   Close the graphic device and clear memory
    #   ===========================================
    dev.off();
    message("R Circos Tick Demo Done ...\n");


