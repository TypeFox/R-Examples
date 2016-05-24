
#       Source code for RCircos package
#   
#       Date created:   January 2, 2013
#       Version:        RCircos.1.0.0
#
#       Last revised on December 1, 2015
#       Version:    RCircos.1.1.3
# 
#       by Hongen Zhang, Ph.D. (hzhang@mail.nih.gov)
#  
#       Genetics Branch
#       Center for Cancer Research 
#       National Cancer Institute
#       National Institutes of Health
#       Bethesda, Maryland 20892
#
# ________________________________________________________________________
# <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>




#   Working environment for RCircos to hold RCircos core components.
#   Last revised on May 16, 2013
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>  

RCircos.Env <- new.env();


#   1.  Print out workflow as a quick guide.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Argument:       None
#   Return value:   None
#
#   Example:    library(RCircos);
#               RCircos.Workflow();
#
RCircos.Workflow <- function()
{
    message("1. Load RCircos library:\n");
    message("   library(RCircos);\n");

    message("2. Load chromosome cytoband data:");
    message("   data(UCSC.HG19.Human.CytoBandIdeogram);");
    message("   cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;\n");

    message("   Other chromosome ideogram data installed:\n");
    message("   UCSC.Mouse.GRCm38.CytoBandIdeogram");
    message("   UCSC.Baylor.3.4.Rat.cytoBandIdeogram\n");

    message("3. Setup RCircos core components:\n");
    message("   RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL, ",
                "10, 0);\n\n");

    message("4. Load input data:\n");
    message("   heatmap.data <- read.table(\"/path/Heatmap.data.txt\", 
                                    sep=\"\\t\", quote=\"\", head=T);");
    message("   hist.data <- read.table(\"/path/histgram.data.txt\", 
                                    sep=\"\\t\", quote=\"\", head=T);");
    message("   link.data <- read.table(\"/path/link.data.txt\", 
                                    sep=\"\\t\", quote=\"\", head=T);\n");

    message("5. Modify plot parameters if necessary:\n");
    message("   rcircos.params <- RCircos.Get.Plot.Parameters()");
    message("   rcircos.params$radiu.len <- 1.5;");
    message("   RCircos.Reset.Plot.Parameters(rcircos.params);\n");

    message("6. Open graphic device:\n");
    message("   RCircos.Set.Plot.Area();\n"); 
    message("   or submit your own code. For example:\n");
    message("   par(mai=c(0.25, 0.25, 0.25, 0.25));");
    message("   plot.new();");
    message("   plot.window(c(-2.5,2.5), c(-2.5, 2.5));\n");

    message("7. Call plot function to plot each data track:\n");
    message("   RCircos.Chromosome.Ideogram.Plot();");
    message("   RCircos.Heatmap.Plot(heatmap.data, data.col=5, 
                track.num=1, side=\"in\");");
    message("   RCircos.Histogram.Plot(hist.data, data.col=4, 
                track.num=4, side=\"in\");");
    message("   RCircos.Link.Plot(link.data, track.num=5, 
                by.chromosome=FALSE);\n");

    message("8. Close the graphic device if you was plotting to file:\n");
    message("   dev.off();\n");
}


#    2. Initialize RCircos core components
#   _________________________________________________________________________
#    <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#  
#    Argument:    
#
#       cyto.info
#       A data frame which contains chromosome ideogram data, e.g., an object
#       returned from calls to the function of read.table() which read a file 
#       containing all information in the cytoBandIdeo table querried from 
#       UCSC genome browser.
#
#       chr.exclude
#       Character vector, chromosome name(s) to be excluded, e.g., "chrM"
#
#       tracks.inside and tracks.outside
#       Non-negative integer, number of data tracks inside or outside of 
#       chromosome ideogram
#
#    Return value:    None
#
#    Example:   library(RCircos);
#               RCircos.Set.Core.Components(cyto.info, 
#                   chr.exclude=NULL, 10, 0);
#

RCircos.Set.Core.Components <- function(cyto.info, chr.exclude=NULL, 
                tracks.inside, tracks.outside)
{
    if(tracks.inside < 0 || tracks.outside < 0) 
    { stop("Track number cannot be smaller than 0.\n"); }

    #   Step 1. validate cyto.info for correct chromosome start and
    #   end positions of each chromosome band. The data will not be 
    #   held for the RCircos environment
    #   ================================================================
    # 
    cyto.band.data <- RCircos.Validate.Cyto.Info(cyto.info, chr.exclude);

    #    Step 2. Initialize RCircos core components
    #   ========================================================
    #
    RCircos.Initialize.Parameters(tracks.inside,tracks.outside);
    RCircos.Set.Cytoband.data(cyto.band.data);
    RCircos.Set.Base.Plot.Positions();
}


#
#   3.  Get methods for retrieving RCircos core components stored in RCircos 
#       environment. These function are used for both outside or internally
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Argument:    None
#
#   Return value:   Paramters, cytoband, or positions stored 
#                   in RCircos environment.
#
#   Example:    plot.param <- RCircos.Get.Plot.Parameters();
#               plot.cyto  <- RCircos.Get.Plot.Ideogram();
#               plot.pos   <- RCircos.Get.Plot.Positions();
#
RCircos.Get.Plot.Parameters <- function()
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    return (RCircosEnvironment[["RCircos.PlotPar"]]);
}

RCircos.Get.Plot.Ideogram <- function()
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    return (RCircosEnvironment[["RCircos.Cytoband"]]);
}

RCircos.Get.Plot.Positions <- function()
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    return (RCircosEnvironment[["RCircos.Base.Position"]]);
}


#   4.  Print out all paramters. This could be ran anytime for checking the 
#       current values of parameters
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Argument:       None
#   Return value:   None
#
#   Example:    RCircos.List.Parameters();
#
RCircos.List.Parameters <- function()
{
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    message("Parameters for current RCircos session.\n\n");

    #   parameters relative to radius.len
    #   ===============================================================
    #
    message("Parameters relative to radius.len:");
    message("==================================\n");

    message(paste("radius.len:\t",     RCircos.Par$radius.len));
    message(paste("chr.ideog.pos:\t",  RCircos.Par$chr.ideog.pos));
    message(paste("highlight.pos:\t",  RCircos.Par$highlight.pos));
    message(paste("chr.name.pos:\t",   RCircos.Par$chr.name.pos));
    message(paste("plot.radius:\t",    RCircos.Par$plot.radius, "\n"));

    message(paste("track.in.start:\t", RCircos.Par$track.in.start));
    message(paste("track.out.start:",  RCircos.Par$track.out.start, "\n"));

    message(paste("chrom.width:\t",    RCircos.Par$chrom.width));
    message(paste("track.padding:\t",  RCircos.Par$track.padding));
    message(paste("track.height:\t",   RCircos.Par$track.height, "\n"));

    #   parameters relative to chromosome unit
    #   =================================================================
    #
    message("Parameters relative to chromosome unit:");
    message("=======================================\n");
    
    message(paste("base.per.unit:\t",  RCircos.Par$base.per.unit));
    message(paste("chrom.paddings:\t", RCircos.Par$chrom.paddings));
    message(paste("heatmap.width:\t",  RCircos.Par$heatmap.width));
    message(paste("hist.width:\t",     RCircos.Par$hist.width, "\n"));

    #   General R graphic parameters 
    #   =================================================================
    #
    message("General R graphic parameters:");
    message("=============================\n");
    
    message(paste("text.size:\t",      RCircos.Par$text.size));
    message(paste("highlight.width:",  RCircos.Par$highlight.width));
    message(paste("point.type:\t",     RCircos.Par$point.type));
    message(paste("point.size:\t",     RCircos.Par$point.size, "\n"));

    message(paste("text.color\t",       RCircos.Par$text.color));
    message(paste("heatmap.color\t",   RCircos.Par$heatmap.color));
    message(paste("hist.color:\t",     RCircos.Par$hist.color));
    message(paste("line.color:\t",     RCircos.Par$line.color));
    message(paste("scatter.color:\t",  RCircos.Par$scatter.color));
    message(paste("tile.color:\t",     RCircos.Par$tile.color, "\n"));

    message(paste("track.background:\t", RCircos.Par$track.background));
    message(paste("grid.line.color:\t", RCircos.Par$grid.line.color, "\n"));

    message(paste("Bezier.point:\t",   RCircos.Par$Bezier.point));
    message(paste("max.layers:\t",     RCircos.Par$max.layers))
    message(paste("sub.tracks:\t",     RCircos.Par$sub.tracks, "\n"))

    #    User friendly note
    #   =================================================================
    #
    message("Following are an example to change RCircos plot parameters:");
    message("===========================================================\n");

    message("params <- RCircos.Get.Plot.Parameters();");  
    message("params$radius.len <- 2.0;");
    message("params$base.per.unit <- 5000;");
    message("RCircos.Reset.Plot.Parameters(params)\n");

    message("Chromosome ideogram data were automatically modified.\n");
}


#   5.  Open a new device with the build-in plot.radius. If user wants an  
#       image file, correct image type must be created and closed from 
#       command line or other script. This function can also be replaced 
#       with a script from command line if user known how much of the plot 
#       area are needed.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#   
#    Argument:      None
#    Return value:  None
#
#    Example:    RCircos.Set.Plot.Area();
#

RCircos.Set.Plot.Area <- function()
{
    RCircos.Par <- RCircos.Get.Plot.Parameters()

    par(mai=c(0.25, 0.25, 0.25, 0.25));
    plot.new();
    plot.window(c(-1*RCircos.Par$plot.radius, RCircos.Par$plot.radius), 
        c(-1*RCircos.Par$plot.radius, RCircos.Par$plot.radius));
}


#   6.  Reset RCircos plot parameters
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Argument:
#   
#   new.params: 
#       A set of new plot parameters came from outside and is not the one
#       stored in RCircos.Env environment.
#
#   Return value:   None
#
#   Example:    plot.params <- RCircos.Get.Plot.Parameters();
#               plot.params$radius.len <- 2;
#               RCircos.Reset.Plot.Parameters(plot.params);
#
RCircos.Reset.Plot.Parameters <- function(new.params)
{
    #   Check parameters for negative values 
    #   =================================================
    #
    point.chr       <- new.params$point.type;
    text.color      <- new.params$text.color;
    heatmap.color   <- new.params$heatmap.color;
    hist.color      <- new.params$hist.color;
    line.color      <- new.params$line.color;
    scatter.color   <- new.params$scatter.color;
    tile.color      <- new.params$tile.color;
    bg.color        <- new.params$track.background;
    grid.color      <- new.params$grid.line.color;

    params <- unlist(new.params);
    params <- params[-which(params == point.chr)];
    params <- params[-which(params == text.color)];
    params <- params[-which(params == heatmap.color)];
    params <- params[-which(params == hist.color)];
    params <- params[-which(params == line.color)];
    params <- params[-which(params == scatter.color)];
    params <- params[-which(params == tile.color)];
    params <- params[-which(params == bg.color)];
    params <- params[-which(params == grid.color)];

    params <- as.numeric(params);
    if(sum(is.na(params)) > 0) 
    { stop("Plot parameters except of point.type must be numeric.\n"); }

    if(sum(params < 0) > 0) 
    { stop("Plot parameters cannot have negative values.\n"); }

    #    Hold old params for later use
    #   =================================================
    #
    old.params <- RCircos.Get.Plot.Parameters();

    #   Remove the old location data from cyto.band.data 
    #   since it will be recalculated and added.
    #   =================================================

    cyto.band.data <- RCircos.Get.Plot.Ideogram();
    cyto.band.data$Unit <- NULL;
    cyto.band.data$Location <- NULL;    

    #   Update chromosome padding based on number of base.per.unit.
    #   Padding length should be always no more than 3/1000 total 
    #   chromosome units
    #   ===========================================================
    #
    if(new.params$chrom.paddings!=0) 
    { 
        padding.const <- 3000/1000000;
        genome.lenth <- sum(cyto.band.data$Length);
        total.units <- genome.lenth/new.params$base.per.unit;
        the.padding <- round(padding.const*total.units, digits=0);

        if(new.params$chrom.paddings>the.padding) {
            message(paste("\nNote: chrom.padding", 
                new.params$chrom.paddings, "is too big,",
                "and was reset to", the.padding, "\n"));
            new.params$chrom.paddings <- the.padding; 
        }
    }

    #    Always keep radius.len 1.0 or greater.
    #   =================================================

    if(new.params$radius.len < 1.0)
    { 
        message("Note: radius.len needs be at least 1.0\n\n"); 
        new.params$radius.len <- 1.0;
    }

    #    Parameters for ideogram are binded to radius.len:
    #    chr.ideog.pos, highlight.pos, chr.name.pos, 
    #    track.in.start, track.out.start, and plot.radius
    #   =================================================
    #
    new.params$chr.ideog.pos <- new.params$radius.len + 0.1;
    new.params$highlight.pos <-new.params$chr.ideog.pos + 
                        new.params$chrom.width + 0.05;
    new.params$chr.name.pos  <- new.params$highlight.pos + 0.15;
    new.params$highlight.width=round(new.params$radius.len, digits=0);

    new.params$track.in.start <- new.params$radius.len; 
    new.params$track.out.start <- new.params$chr.name.pos+0.2;
    
    differ <- old.params$plot.radius - old.params$radius.len;
    new.params$plot.radius <- new.params$radius.len + differ;

    #    Remove all objects from RCircos environment
    #   =======================================================
    #
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());

    RCircosEnvironment[["RCircos.PlotPar"]] <- NULL;
    RCircosEnvironment[["RCircos.Cytoband"]] <- NULL;
    RCircosEnvironment[["RCircos.Base.Position"]] <- NULL;

    #   Put the plot parameter in RCircos environment and reset
    #   cytoband data and base positions
    #   =================================================
    #
    RCircosEnvironment[["RCircos.PlotPar"]] <- new.params;
    RCircos.Set.Cytoband.data(cyto.band.data);
    RCircos.Set.Base.Plot.Positions();
}


#   7.  Draw chromosome ideogram. Graphic device must be initialized first
#   ______________________________________________________________________
#    <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#    Argument:      None
#    Return value:  None
#
#    Example:    RCircos.Chromosome.Ideogram.Plot();
#

RCircos.Chromosome.Ideogram.Plot <- function()
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Pos  <- RCircos.Get.Plot.Positions();
    RCircos.Par  <- RCircos.Get.Plot.Parameters();

    #    Plot chromosome outlines, chromosome names, and 
    #    chromosome highlights
    #   =================================================
    #
    outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width;
    inner.location <- RCircos.Par$chr.ideog.pos;

    chroms <- unique(RCircos.Cyto$Chromosome);
    for(a.chr in 1:length(chroms))
    {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
        start   <- the.chr$Location[1]- the.chr$Unit[1] + 1;
        end     <- the.chr$Location[nrow(the.chr)];
        mid     <- round((end-start+1)/2, digits=0)+start;

        chr.color <- the.chr$ChrColor[nrow(the.chr)];

        #    Draw chromosome outlines
        #   =================================================
        #
        pos.x<- c(RCircos.Pos[start:end,1]*outer.location, 
                  RCircos.Pos[end:start,1]*inner.location);
        pos.y<- c(RCircos.Pos[start:end,2]*outer.location, 
                RCircos.Pos[end:start,2]*inner.location);
        polygon(pos.x, pos.y);

        #    Add chromosome names
        #   =================================================
        #
        chr.name <- sub(pattern="chr", replacement="", chroms[a.chr]);
        text(RCircos.Pos[mid,1]*RCircos.Par$chr.name.pos,
            RCircos.Pos[mid,2]*RCircos.Par$chr.name.pos,
            label=chr.name, srt=RCircos.Pos$degree[mid]);

        #    Add chromosome highlights
        #   =================================================
        #
        lines(RCircos.Pos[start:end,]*RCircos.Par$highlight.pos, 
                col=chr.color, lwd=RCircos.Par$highlight.width);
    }

    #    Add chromosome bands (Giema stain positive only)
    #   =================================================
    #
    for(a.band in 1:nrow(RCircos.Cyto))
    {
        a.color <- RCircos.Cyto$BandColor[a.band];
        if(a.color=="white") { next; }

        start <- RCircos.Cyto$Location[a.band] - RCircos.Cyto$Unit[a.band] + 1;
        end   <- RCircos.Cyto$Location[a.band];

        pos.x<- c(RCircos.Pos[start:end,1]*outer.location, 
                RCircos.Pos[end:start,1]*inner.location);
        pos.y<- c(RCircos.Pos[start:end,2]*outer.location, 
                RCircos.Pos[end:start,2]*inner.location);
        polygon(pos.x, pos.y, col=a.color, border=NA);
    }
}


#   8.  Draw connectors between chromosome ideogram and gene lables. The  
#       argument (connect.data) has two columns for paired point locations  
#       on two tracks(gene label and chromosomes ideogram). The first 
#       column is outer points, and the second column is inner points. 
#
#       The points are sorted by relative positions on their chromosomes  
#       and are held in the last two columns of connect.data
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       genomic.data
#       A data frame with the first four columns for chromosome name, start 
#       position, end position, and names of genes.
#
#       track.num
#       Non-negative integer, the track number from chromosome ideogram.
#
#       side
#       Character vector, must be either "in" or "out".
#
#    Return value:  None
#    Example:       RCircos.Gene.Connector.Plot(genomic.data, 1, "in")
#

RCircos.Gene.Connector.Plot <- function(genomic.data, track.num, side)
{
    if(track.num<1) { stop("Track number cannot be smaller than 1.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Construct Connector data from gene name data
    #   =========================================================
    #
    gene.data    <- RCircos.Get.Plot.Data(genomic.data, "plot");
    label.data   <- RCircos.Get.Gene.Label.Locations(gene.data);
    connect.data <- data.frame(label.data$Location, label.data$Label.Position);

    #   Plot position for current track. 
    #   =========================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);
    out.pos   <- locations[1];
    in.pos    <- locations[2];

    #   Heights for the two vertical lines of connectors and
    #   the horizontal line range
    #   =========================================================
    #
    v.height <- round((out.pos-in.pos)/10, digits=4);
    h.range <- out.pos - in.pos - 2*v.height;

    top.loc <- out.pos - v.height;
    bot.loc <-  in.pos +  v.height;

    #   Connector colors
    #   =========================================================
    #
    line.colors <- RCircos.Get.Plot.Colors(label.data, RCircos.Par$text.color);

    #   To handle connector plot on both sides
    #   =========================================================
    #
    if(side=="in") {
        genomic.col <- ncol(connect.data) - 1;
        label.col <- ncol(connect.data);
    } else {
        genomic.col <- ncol(connect.data);
        label.col <- ncol(connect.data) - 1;
    }

    chroms <- unique(connect.data[,1]);
    for(a.chr in 1:length(chroms))
    {
        chr.row <- which(connect.data[,1]==chroms[a.chr]);
        total <- length(chr.row);

        for(a.point in 1:total)
        {
            p1 <- connect.data[chr.row[a.point], genomic.col];
            p2 <- connect.data[chr.row[a.point], label.col];

            #   draw top vertical line
            #   ============================================================
            #
            lines(c(RCircos.Pos[p1, 1]*out.pos, RCircos.Pos[p1, 1]*top.loc),
                  c(RCircos.Pos[p1, 2]*out.pos, RCircos.Pos[p1, 2]*top.loc),
                  col=line.colors[chr.row[a.point]]);

            #   draw bottom vertical line
            #   ===========================================================
            #
            lines(c(RCircos.Pos[p2, 1]*bot.loc, RCircos.Pos[p2, 1]*in.pos),
                  c(RCircos.Pos[p2, 2]*bot.loc, RCircos.Pos[p2, 2]*in.pos),
                  col=line.colors[chr.row[a.point]]);

            #   draw horizontal line
            #   ============================================================
            #
            lines(c(RCircos.Pos[p1, 1]*top.loc, RCircos.Pos[p2, 1]*bot.loc),
                  c(RCircos.Pos[p1, 2]*top.loc, RCircos.Pos[p2, 2]*bot.loc),
                  col=line.colors[chr.row[a.point]]);
        }
    }
}


#   9.  Label genes beside of track. This is only suitable for small number 
#       of labels. When cex=0.4, each character of the lable will occupy 
#       about 5000 units. This is the best visulization for a 8x8 inche image.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       gene.data
#       A data frame returned from RCircos.Get.Gene.Label.Locations 
#       (genomic.data). The genomic.data has three leading columns 
#       for genomic positions followed by gene names.
#
#       name.col
#       Non-negative integer, number of column in gene.data for data 
#       (gene names) to be plotted 
#
#       track.num
#       Non-negative integer, number of the track from chromosome ideogram.
#
#       side
#       Character vector, must be either "in" or "out"
#
#   Return value:    None
#
#   Example:    RCircos.Gene.Label(gene.data, 4, 3, "in")
#

RCircos.Gene.Name.Plot <- function(gene.data, name.col, track.num, side)
{
    if(track.num < 1) { stop("Track number cannot be smaller than 0.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("Side must be either in or out.\n"); }

    if(name.col < 4) 
    { stop("Column number for gene names must be 4 or bigger.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Convert raw data to plot data. The raw data will be 
    #   validated first during the convertion
    #   =====================================================
    #
    gene.data <- RCircos.Get.Plot.Data(gene.data, "plot");
    gene.data <- RCircos.Get.Gene.Label.Locations(gene.data);

    #    Label positions
    #   ===================================================
    #
    right.side <- nrow(RCircos.Pos)/2;
    locations <- RCircos.Track.Positions(side, track.num);
    the.points <- as.numeric(gene.data[, ncol(gene.data)]);

    if(side=="in") {
        label.pos <- locations[1]; 
        text.side <- rep(4, nrow(gene.data));
        text.side[the.points <= right.side] <- 2;    
    } else {
        label.pos  <- locations[2];
        text.side <- rep(2, nrow(gene.data));
        text.side[the.points<=right.side] <- 4; 
    }

    #   Label colors
    #   ===================================================
    #
    text.colors <- RCircos.Get.Plot.Colors(gene.data, RCircos.Par$text.color);

    #    Plot labels
    #   ===================================================
    #
    for(a.text in 1:nrow(gene.data))
    {
        gene.name <- as.character(gene.data[a.text, name.col]);
        rotation <- RCircos.Pos$degree[the.points[a.text]];

        text(RCircos.Pos[the.points[a.text], 1]*label.pos,
             RCircos.Pos[the.points[a.text], 2]*label.pos,
            label=gene.name, pos=text.side[a.text], cex=RCircos.Par$text.size, 
            srt=rotation, offset=0, col=text.colors[a.text]);
    }
}


#   10. Draw one track of heatmap with blue and read colors. The first 
#       four columns of headmap data must be chromosome, chromStart, 
#       chromEnd, and gene names.
#   _________________________________________________________________________
#    <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
# 
#    Arguments:
#
#       heatmap.data
#       A data frame with returned from RCircos.Get.Plot.Data(genomic.data, 
#       plot.type). Heatmap data must have four leading columns for gene 
#       position and gene names.
#
#       data.col
#       Non-negative integer, number of column in heatmap.data for data to 
#       be plotted.
#
#       track.num
#       Non-negative integer, number of the track from chromosome ideogram.
#
#       side
#       Character vector, must be either "in" or "out".
#
#    Return value:  None
#
#    Example:   RCircos.Heatmap.Plot(heatmap.data, 3, 3, "in")
#
RCircos.Heatmap.Plot <- function(heatmap.data, data.col, track.num, side)
{
    if(track.num<1) { stop("Track number cannot be smaller than 0.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    if(data.col<5) { stop("Heatmap data column must be 5 or bigger.\n"); }

    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Convert raw data to plot data. The raw data will be validated
    #   first during the convertion
    #   =============================================================
    #
    heatmap.data <- RCircos.Get.Plot.Data(heatmap.data, "plot");

    #   Colors for different data values
    #   =============================================================
    #
    ColorRamp <- RCircos.Get.Heatmap.ColorScales(RCircos.Par$heatmap.color);

    #   Color level. Heatmap data has to have four leading columns
    #   for genomic position and gene(lable) names. Also color 
    #   level must be calculated with all data columns.
    #   =============================================================
    #
    columns    <- 5:(ncol(heatmap.data)-1);
    min.value  <- min(as.matrix(heatmap.data[, columns]));
    max.value  <- max(as.matrix(heatmap.data[, columns]));
    ColorLevel <- seq(min.value, max.value, length=length(ColorRamp));

    #    Each heatmap cell is centered on data point location. Make
    #    sure each one will be in the range of it chromosome
    #   =============================================================
    #
    heatmap.locations <- as.numeric(heatmap.data[, ncol(heatmap.data)]);
    start <- heatmap.locations - RCircos.Par$heatmap.width/2;
    end   <- heatmap.locations + RCircos.Par$heatmap.width/2;

    data.chroms <- as.character(heatmap.data[,1]);
    chromosomes <- unique(data.chroms);
    cyto.chroms <- as.character(RCircos.Cyto$Chromosome);

    for(a.chr in 1:length(chromosomes))
    {
        cyto.rows <- which(cyto.chroms==chromosomes[a.chr]);
        locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]);
        chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]];
        chr.end   <- max(locations);

        data.rows <- which(data.chroms==chromosomes[a.chr]);
        start[data.rows[start[data.rows]<chr.start]] <- chr.start;
        end[data.rows[end[data.rows]>chr.end]] <- chr.end;
    }

    #    Plot position for current track. 
    #   ==================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);
    out.pos <- locations[1];
    in.pos  <- locations[2];

    #    outline of chromosomes. No lines inside.
    #   =========================================
    #
    chroms <- unique(RCircos.Cyto$Chromosome);
    for(a.chr in 1:length(chroms))
    {
        the.chr   <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr], ];
        the.start <- the.chr$Location[1] - the.chr$Unit[1] + 1;
        the.end   <- the.chr$Location[nrow(the.chr)];

        polygon.x<- c(RCircos.Pos[the.start:the.end, 1]*out.pos, 
                      RCircos.Pos[the.end:the.start, 1]*in.pos);
        polygon.y<- c(RCircos.Pos[the.start:the.end, 2]*out.pos, 
                      RCircos.Pos[the.end:the.start, 2]*in.pos);
        polygon(polygon.x, polygon.y, col="white");
    }

    #    Plot heatmap for each gene.
    #   =================================================
    #
    heatmap.value <- as.numeric(heatmap.data[, data.col]);
    for(a.point in 1:length(heatmap.value))
    {
        the.level <- which(ColorLevel>=heatmap.value[a.point]);
        cell.color <- ColorRamp[min(the.level)];
        
        the.start <- start[a.point];
        the.end <- end[a.point];

        polygon.x<- c(RCircos.Pos[the.start:the.end,1]*out.pos, 
                RCircos.Pos[the.end:the.start,1]*in.pos);
        polygon.y<- c(RCircos.Pos[the.start:the.end,2]*out.pos, 
                RCircos.Pos[the.end:the.start,2]*in.pos);
        polygon(polygon.x, polygon.y, col=cell.color, border=NA);
    }
}


#   11.    Draw one track of histogram. 
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       hist.data
#       A data frame with returned from RCircos.Get.Plot.Data(genomic.data, 
#       plot.type) which has three leading columns for genomic positions.
#
#       data.col
#       Non-negatuve integer, number of column in heatmap.data for data to 
#       be plotted
#
#       track.num
#       Non-negatuve integer, number of the track from chromosome ideogram.
#
#       side
#       Character vector, must be either "in" or "out".
#
#    Return value:  None
#    Example:       RCircos.Histogram.Plot(hist.data, 4, 2, "in")
#
RCircos.Histogram.Plot <- function(hist.data, data.col, track.num, side)
{
    #   Arguments validate
    #   =============================================================
    #
    if(track.num<1) { stop("Track number cannot be smaller than 1.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    if(data.col<4) { stop("Heatmap data column must be 4 or bigger.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    #    Convert raw data to plot data. The raw data will be validated
    #    first during the convertion
    #   =============================================================
    #
    hist.data <- RCircos.Get.Plot.Data(hist.data, "plot");

    #    Each Histogram cell is centered on data point location. Make
    #    sure each one will be in the range of it chromosome
    #   =============================================================
    #
    hist.locations <- as.numeric(hist.data[, ncol(hist.data)]);
    start <- hist.locations - RCircos.Par$hist.width;
    end   <- hist.locations + RCircos.Par$hist.width;

    data.chroms <- as.character(hist.data[,1]);
    chromosomes <- unique(data.chroms);
    cyto.chroms <- as.character(RCircos.Cyto$Chromosome);

    for(a.chr in 1:length(chromosomes))
    {
        cyto.rows <- which(cyto.chroms==chromosomes[a.chr]);
        locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]);
        chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]];
        chr.end   <- max(locations);

        data.rows <- which(data.chroms==chromosomes[a.chr]);
        start[data.rows[start[data.rows]<chr.start]] <- chr.start;
        end[data.rows[end[data.rows]>chr.end]] <- chr.end;
    }

    #    Plot position for current track. 
    #   =============================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);
    out.pos <- locations[1];
    in.pos <- locations[2];

    #       Label colors
    #   =============================================================
    #
    hist.colors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color); 

    #    Draw histogram
    #   =============================================================
    #
    num.subtrack <- RCircos.Par$sub.tracks;
    RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks, NULL);

    for(a.point in 1:nrow(hist.data))
    {
        hist.height <- hist.data[a.point, data.col];
        the.start <- start[a.point];
        the.end <- end[a.point];

        #    Plot rectangle with specific height for each 
        #    data point
        #   =============================================================
        #
        height <- in.pos + RCircos.Par$track.height*hist.height;

        polygon.x<- c(RCircos.Pos[the.start:the.end, 1]*height, 
                      RCircos.Pos[the.end:the.start, 1]*in.pos);
        polygon.y<- c(RCircos.Pos[the.start:the.end, 2]*height, 
                      RCircos.Pos[the.end:the.start, 2]*in.pos);
        polygon(polygon.x, polygon.y, col=hist.colors[a.point], border=NA);
    }
}


#   12.    Draw one track of line plot
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#    Arguments:
#
#       line.data
#       A data frame with returned from RCircos.Get.Plot.Data(genomic.data). 
#       line.data has three leading columns for genomic positions.
#
#       data.col
#       Non-negative integer, number of column in heatmap.data for data to 
#       be plotted.
#
#       track.num
#       Non-negative integer, number of the track from chromosome ideogram
#
#       side
#       Character vector, must be either "in" or "out"
#
#    Return value:  None
#    Example:       RCircos.Line.Plot(line.data, 4, 3, "in")
#
RCircos.Line.Plot <- function(line.data, data.col, track.num, side)
{
    #   Arguments validate
    #   =========================================================
    #
    if(track.num<1) 
    { stop("Track number cannot be smaller than 1.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    if(data.col<4) 
    { stop("heatmap data column must be 4 or bigger.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Convert raw data to plot data. The raw data will be 
    #   validated first during the convertion
    #   ===================================================
    #
    line.data <- RCircos.Get.Plot.Data(line.data, "plot");

    #    Plot position for current track. 
    #   ===================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);
    out.pos <- locations[1];
    in.pos  <- locations[2];

    #    Check if the data has negative values such as copy number
    #    change or log2 of fold change. If yes the zero line will use 
    #    the middle of track height otherwise the inner boundary
    #   ==========================================================
    #
    if(min(as.numeric(line.data[,data.col]))>=0) {
        point.bottom <- in.pos; data.ceiling <- 10;
    } else {  
        point.bottom <- in.pos + (RCircos.Par$track.height/2);  
        data.ceiling <- 5;
    }
    sub.height <- out.pos-point.bottom;

    #       Line colors
    #   ===================================================
    #
    line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color); 

    #    Start plotting. Line plot is connecting two neighbor points
    #    so no exception catch needed.
    #   ===================================================
    #
    RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks, NULL);

    for(a.point in 1:(nrow(line.data)-1))
    {
        point.one <-  line.data[a.point, ncol(line.data)];
        point.two <- line.data[a.point+1, ncol(line.data)];

        #    Next data point is on next chromosomes or null. Stop.
        #   ===================================================
        #
        if(line.data[a.point,1]!= line.data[a.point+1,1]) { next;}

        #    Check data value for correct plot ranges
        #   ===================================================
        #
        if(line.data[a.point, data.col] >data.ceiling) { 
            value.one<- data.ceiling;
        } else if(line.data[a.point, data.col] <(-1*data.ceiling)) {  
            value.one <- data.ceiling*-1;
        } else { value.one <- line.data[a.point, data.col];  }

        height.one <- point.bottom  + value.one/data.ceiling*sub.height;

        if(line.data[a.point+1, data.col] >data.ceiling) { 
            value.two<- data.ceiling;
        } else if(line.data[a.point+1, data.col] <(-1*data.ceiling)) {  
            value.two <- data.ceiling*-1;
        } else { value.two <- line.data[a.point+1, data.col];  }

        height.two <- point.bottom  + value.two/data.ceiling*sub.height;

        #    Draw lines
        #   ===================================================
        #
        lines(c(RCircos.Pos[point.one, 1]*height.one,
                RCircos.Pos[point.two, 1]*height.two),
              c(RCircos.Pos[point.one, 2]*height.one,
                RCircos.Pos[point.two, 2]*height.two),
                col=line.colors[a.point]);
    }
}


#   13.    Draw one track of scatterplot.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R> 
#
#    Arguments:
#
#       scatter.data
#       A data frame returned from RCircos.Get.Plot.Data(genomic.data, 
#       plot.type). The scatter.data has three leading columns for 
#       genomic positions.
#
#       data.col
#       Non-negative integer, number of column in heatmap.data for data 
#       to be plotted
#
#       track.num
#       Non-negative integer, number of the track from chromosome ideogram.
#
#       side
#       Character vector, must be either "in" or "out"
#
#       by.fold
#       Positive integer. Use for color cotrol of the points.
#
#   Return value:   None
#   Example:        RCircos.Scatter.Plot(scatter.data, 5, 3, "in", 1)
#
RCircos.Scatter.Plot <- function(scatter.data, data.col, track.num, side, 
        by.fold=0)
{
    #   Arguments validate
    #   =================================================================
    #
    if(track.num<1) { stop("Track number cannot be smaller than 1.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    if(data.col<4) { stop("Heatmap data column must be 4 or bigger.\n"); }
    if(by.fold<0)  { stop("Value of by.fold cannot be negative"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #    Convert raw data to plot data. The raw data will be 
    #    validated first during the convertion
    #   ========================================================
    #
    scatter.data <- RCircos.Get.Plot.Data(scatter.data, "plot");

    #    Plot position for current track. 
    #   ========================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);
    out.pos <- locations[1];
    in.pos  <- locations[2];

    #    Check if the data has negative value such as copy number
    #    change or log2 of fold change. If yes the zero line will use 
    #    the middle of track height otherwise the inner boundary
    #   ========================================================
    #
    if(min(as.numeric(scatter.data[,data.col]))>=0) {
        point.bottom <- in.pos; data.ceiling <- 10;
    } else {  
        point.bottom <- in.pos + (RCircos.Par$track.height/2);  
        data.ceiling <- 5;
    }
    sub.height <- out.pos-point.bottom;

    #       scatter point colors
    #   ========================================================
    #
    scatter.colors <- RCircos.Get.Plot.Colors(scatter.data, 
                            RCircos.Par$scatter.color); 

    #    Start plotting
    #   ========================================================
    #
    RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks, NULL);

    for(a.point in 1:nrow(scatter.data))
    {
        the.point <- scatter.data[a.point, ncol(scatter.data)];

        #    Adjust the data value to avoid overflow
        #   ========================================================
        #
        if(scatter.data[a.point, data.col] >data.ceiling) { 
            the.value<- data.ceiling;
        } else if(scatter.data[a.point, data.col] <(-1*data.ceiling)) {  
            the.value <- data.ceiling*-1;
        } else { the.value <- scatter.data[a.point, data.col]; }

        #    Color for the point
        #   ========================================================
        #
        color <- scatter.colors[a.point];
        if(by.fold>0) {
            if(the.value>=by.fold) { color <- "red"; 
            } else if (the.value<=-by.fold) { color <- "blue";
            } else {  color <- "black"; }
        } 

        #    plot a scatter
        #   ========================================================
        #
        height <- point.bottom  + the.value/data.ceiling*sub.height;
        points(RCircos.Pos[the.point, 1]*height,
                RCircos.Pos[the.point ,2]*height,
                col=color, pch=RCircos.Par$point.type, 
                cex=RCircos.Par$point.size);
    }
}


#   14. Draws one track of tiles. Note: Tile plot needs genomic position 
#       only and data column is not requied.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#    Arguments:
#
#       tile.data
#       A data frame returned from RCircos.Get.Plot.Data(genomic.data, 
#       cyto.band). The tile.data has three leading columns for genomic 
#       positions.
#
#       track.num
#       Non-negative integer, number of the track from chromosome ideogram.
#
#       side
#       Character vector, must be either "in" or "out"
#
#    Return value:  None
#    Example:       RCircos.Tile.Plot(tile.data, track.num=3, side="in", 1)
#
RCircos.Tile.Plot <- function(tile.data, track.num, side)
{
    #   Arguments validate
    #   =================================================
    # 
    if(track.num<1) 
    { stop("Track number cannot be smaller than 1.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #    Convert raw data to plot data. The raw data will be 
    #    validated first during the convertion
    #   ====================================================
    #
    tile.data <- RCircos.Get.Plot.Data(tile.data, "plot");

    #    Assign a layer number to each data point and find 
    #    the maxium layer number
    #   =================================================
    #
    the.layer <- 1;
    the.chr <- tile.data[1,1];
    start <- tile.data[1,2];  
    end  <- tile.data[1,3];
    tile.layers <- rep(1, nrow(tile.data));

    for(a.row in 2:nrow(tile.data))
    {
        #    Meet a new region without overlap with previous or
        #    a different chromosome, reset relevant variables
        #   =================================================
        #
        if (tile.data[a.row, 2] >= end ) {  
            the.layer <- 1;  
            start <- tile.data[a.row, 2];   
            end <- tile.data[a.row, 3];
        } else if (tile.data[a.row, 1] != the.chr) { 
            the.layer <- 1;  
            the.chr <- tile.data[a.row,1];
            start <- tile.data[a.row, 2]; 
            end <- tile.data[a.row, 3];
        } else {  
            the.layer <- the.layer+1; 
            if(tile.data[a.row, 3]>end) 
            { end <- tile.data[a.row, 3]; }
        }
        tile.layers[a.row] <- the.layer;
    }

    #    Plot position for current track. Adjust them if 
    #    total layers is great than RCircos.Par$max.layers.
    #   =================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);
    out.pos <- locations[1];
    in.pos <- locations[2];

    layer.height <- RCircos.Par$track.height/RCircos.Par$max.layers;
    num.layers <- max(tile.layers);

    if(num.layers>RCircos.Par$max.layers) 
    { 
        if(side=="in")
        {  in.pos <- out.pos - layer.height*num.layers;        
        } else { out.pos <- in.pos + layer.height*num.layers; }

        message(paste("Tiles plot will use more than one track.",
            "Please select correct area for next track.\n"));
    }
    
    if(num.layers<RCircos.Par$max.layers) 
    { layer.height <- RCircos.Par$track.height/num.layers; }

    #       Tile colors
    #   =================================================
    #
    tile.colors <- RCircos.Get.Plot.Colors(tile.data, RCircos.Par$tile.color); 

    #    Start plotting
    #   =====================================================
    #
    RCircos.Track.Outline(out.pos, in.pos, num.layers, NULL);

    the.loc <- ncol(tile.data);
    for(a.row in 1:nrow(tile.data))
    {
        tile.len <- tile.data[a.row,3] - tile.data[a.row,2];
        tile.range <- round(tile.len/RCircos.Par$base.per.unit/2, digits=0);
        start <- tile.data[a.row,the.loc] - tile.range
        end   <- tile.data[a.row,the.loc] + tile.range;

        layer.bot <- in.pos + layer.height*(tile.layers[a.row]-1);
        layer.top <- layer.bot + layer.height*0.8;

        polygon.x<- c(RCircos.Pos[start:end,1]*layer.top, 
                RCircos.Pos[end:start,1]*layer.bot);
        polygon.y<- c(RCircos.Pos[start:end,2]*layer.top, 
                RCircos.Pos[end:start,2]*layer.bot);
        polygon(polygon.x, polygon.y, col=tile.colors[a.row]);
    }
}


#   15. Draw link lines between two chromosomes. Link is always inside of 
#       ideogram and data does not need to be converted with function 
#       RCircos.Get.Plot.Data(data, plot.type)
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#  
#    Arguments:
#
#       link.data
#       A data frame with paired genomic positions in each row.
#
#       track.num
#       Non-negative integer, number of the track from chromosome ideogram.
#
#       by.chromosome
#       Logical, if true, intrachromosome will be red, otherwise random 
#       rainbow color will be used
#
#   Return value:   None
#   Example:        RCircos.Link.Plot(link.data, 9, FALSE);
#
RCircos.Link.Plot <- function(link.data, track.num, by.chromosome=FALSE)
{
    #   Arguments validate
    #   =================================================================
    #
    if(track.num<1) { stop("Track number cannot be smaller than 1.\n"); }
    if(by.chromosome!=TRUE && by.chromosome!=FALSE)
    { stop("by.chromosome must be either TRUE or FALSE.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Check chromosome names, Start, and End positions
    #   =================================================================
    #
    link.data <- RCircos.Validate.Genomic.Data(link.data, plot.type="link");

    #   Plot position for link track.
    #   =================================================================
    #
    one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;
    start <- RCircos.Par$track.in.start - (track.num-1)*one.track;
    base.positions <- RCircos.Pos*start;

    data.points <- matrix(rep(0, nrow(link.data)*2), ncol=2);
    for(a.link in 1:nrow(link.data))
    {
        data.points[a.link, 1] <- RCircos.Data.Point(
            link.data[a.link, 1], link.data[a.link, 2]);
        data.points[a.link, 2] <- RCircos.Data.Point(
            link.data[a.link, 4], link.data[a.link, 5]);

        if(data.points[a.link, 1]==0 || data.points[a.link, 2]==0)
        {  print("Error in chromosome locations ...");  break; }
    }

    #   Get link line colors for each pair of locations
    #   ============================================================
    #
    link.colors <- RCircos.Get.Link.Colors(link.data, by.chromosome);

    #    Draw link lines for each pair of locations
    #   ===========================================
    #
    for(a.link in 1:nrow(data.points))
    {  
        point.one <- data.points[a.link, 1];
        point.two <- data.points[a.link, 2];
        if(point.one > point.two)
        { 
            point.one <- data.points[a.link, 2];
            point.two <- data.points[a.link, 1];
        }

        P0 <- as.numeric(base.positions[point.one,]);
        P2 <- as.numeric(base.positions[point.two,]);
        links <- RCircos.Link.Line(P0, P2); 
        lines(links$pos.x, links$pos.y, type="l", 
            col=link.colors[a.link] ); 
    }
}


#   16. Validate chromosome ideogram information for correct chromosome order 
#       including of correct order of chromosomes, chromosome start, and 
#       chromosome end positions
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#    Argument:    
#
#       cyto.info
#       A data frame which containing the chromosome ideogram data, e.g.,
#       an object returned by calls to the function of read.table() which
#       read a file containing full information of cytoBandIdeo table 
#       querried from UCSC genome browser. This is not A RCircos component.
#
#       chr.exclude
#       Character vector, chromosome name(s) to be excluded from RCircos plot.
#
#    Return value:  data frame with validated chromosome ideogram data  
#    Example:       RCircos.Validate.Cyto.Info(cyto.info, c("chrX", "chrY")))
#
RCircos.Validate.Cyto.Info <- function(cyto.info, chr.exclude=NULL)
{
    
    #   Set standard headers for cyto.info
    #   =============================================================
    #
    if(ncol(cyto.info)<5) { 
        stop(paste("Cytoband data must have columns for chromosome,",
                "chromStart, chromEnd, band, Stain\n",
                "Current cyto.info columns:", ncol(cyto.info)));  }

    colnames(cyto.info)[1:5] <- c("Chromosome", "ChromStart", 
                                    "ChromEnd", "Band", "Stain");

    #    Any ChromEnd must be greater than its ChromStart
    #   =============================================================
    #
    endGTstart <- cyto.info$ChromEnd>cyto.info$ChromStart;
    if(sum(endGTstart)<length(endGTstart)) 
    { stop("One or more ChromEnd position is less than ChromStart."); }

    #    Check if any chromosome names missing "chr" or any chromosome
    #    names with "_"
    #   =============================================================
    #
    cyto.info$Chromosome <- as.character(cyto.info$Chromosome);
    rows <- grep("chr", cyto.info$Chromosome);

    if(length(rows)< length(cyto.info$Chromosome)) 
    { 
        index <- 1:length(cyto.info$Chromosome);
        index <- index[-rows];
        cyto.info$Chromosome[index] <- paste0("chr", 
            cyto.info$Chromosome[index]); 
    }

    rows <- grep("_", cyto.info$Chromosome);
    if(length(rows)>0) { cyto.info <- cyto.info[-rows,]; }

    #     Check any chromosomes need to be excluded from plot
    #   =============================================================
    #
    chromosomes <- unique(cyto.info$Chromosome);

    if(length(chr.exclude)>0) {

        num.chroms <- grep("chr", chr.exclude);
        if(length(num.chroms)!= length(chr.exclude)) 
        { stop("chromosomes in chr.exclude must have prefix of chr"); } 

        ignore <- chromosomes %in% chr.exclude;
        if(sum(ignore)>0){ chromosomes <- chromosomes[ignore==F]; }
    }

    #    To handle chromosome with character only
    #   =============================================================
    #
    chr.names <- chromosomes[grep("[0-9]", chromosomes)];
    if(length(chr.names) > 0) 
    {
        chr.names <- chr.names[order(as.numeric(sub("chr", "", chr.names)))];

        #    If there are sex chromosomes (X, Y, and others)
        #   =============================================================
        #
        if (length(chr.names) < length(chromosomes)) {

            other.chr <- chromosomes[-grep("[0-9]", chromosomes)];

            #    human chromosome order has to be X, Y, and M
            #    Others have not been handled.
            #   ===============================================
            #
            if(length(grep("M", other.chr)) > 0) {
                chrM <- other.chr[grep("M", other.chr)];
                other.chr <- other.chr[-grep("M", other.chr)];
                other.chr <- other.chr[order(other.chr)];
                other.chr <- c(other.chr, chrM);
            } else {
                other.chr <- other.chr[order(other.chr)];
            }
            chr.names <- c(chr.names, other.chr);
        }
        chromosomes <- chr.names;
    } else {
        chromosomes <- chromosomes[order(chromosomes)];
    }

    #    Check start and end pisitions for each chromosome band. The first
    #    start position for each chromosome must be 0 or 1 and all other
    #    start positions must be greater than its previous end position.
    #   =============================================================
    #
    new.cyto.info <- cyto.info[which(cyto.info$Chromosome==chromosomes[1]),];
    for(a.chr in 1:length(chromosomes))
    {
        the.row <- which(cyto.info$Chromosome==chromosomes[a.chr]);
        the.info <- cyto.info[the.row,];
        the.start <- as.numeric(the.info$ChromStart);
        the.end   <- as.numeric(the.info$ChromEnd);
        the.info <- the.info[order(the.start),];

        if(the.start[1]>1) 
        { stop("Cytoband start position cannot be greater than 1."); }

        #   To handle the genomes without band information 
        #   ==============================================
        #
        if(length(the.row)>1) { 
            for(a.row in 2:nrow(the.info)) {
                if(the.start[a.row] < the.end[(a.row-1)]) 
                { 
                    stop(paste("Cytoband start position",
                        "cannot be less than",
                        "previous end position.")); 
                }
            }
        }

        if(a.chr==1) { new.cyto.info <- the.info; 
        } else { new.cyto.info <- rbind(new.cyto.info, the.info); }
    }

    #   Return the validated cyto.info 
    #   ===============================
    # 
    return (new.cyto.info);
}


#   17. Validate input dataset for correct chromosome names, chromosome 
#       start, and chromosome end positions. Chromosome names will be 
#       converted to character vectors if they are factor variables.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       genomic.data
#       A data frame with genomic position data 
#
#       plot.type
#       Character vector either "plot" or "link"
#
#   Return value:    validated input data frame
#
#   Example:    RCircos.Validate.Genomic.Data(genomic.data, "plot")
#
RCircos.Validate.Genomic.Data <- function(genomic.data, 
            plot.type=c("plot", "link"))
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    #   Plot data has only one chromosome column and link data have two
    #   ===============================================================
    #
    plot.type <- tolower(plot.type);
    if(plot.type=="plot") { 
        chrom.col <- 1; 
    } else if(plot.type=="link") { 
        chrom.col <- c(1,4); 
    } else { stop("Plot type must be \"plot\" or \"line\""); }
   
    for (a.col in 1:length(chrom.col))
    {
        the.col <- chrom.col[a.col];

        #   Make sure chromosome names in genomic data have prefix
        #   =========================================================
        # 
        genomic.data[,the.col] <- as.character(genomic.data[,the.col]);
        for(a.row in 1:nrow(genomic.data)) {
            if(length(grep("chr", genomic.data[a.row,the.col])) == 0) 
            { genomic.data[a.row,the.col] <- paste("chr", 
                genomic.data[a.row,the.col], sep=""); }
        }

        #   Make sure chromosomes in input data are all included 
        #   in chromosome ideogram data
        #   =========================================================
        #
        cyto.chroms <- unique(as.character(RCircos.Cyto$Chromosome));
        data.chroms <- unique(as.character(genomic.data[,the.col]));
        if(sum(data.chroms %in% cyto.chroms) < length(data.chroms)) 
        { 
            message(paste("Some chromosomes are in genomic data only",
                            "and have been removed.\n\n"));

            all.chroms <- as.character(genomic.data[,the.col]);
            genomic.data <- genomic.data[all.chroms %in% cyto.chroms,];
        }
        data.chroms <- unique(as.character(genomic.data[,the.col]));

        #   Make sure chromosome start and end postions in genomic
        #   data are not negative.
        #   ========================================================
        # 
        if(min(genomic.data[,the.col+1]) < 0) 
        { stop("Error! chromStart position less than 0."); }
        if(min(genomic.data[,the.col+2]) < 0) 
        { stop("Error! chromEnd position less than 0.");  }    

        #   Make sure chromosome start and end locations in genomic
        #   data are not out of chromosome length
        #   ========================================================
        #
        for(a.chr in 1:length(data.chroms))
        {
            data.chr <- data.chroms[a.chr];
            in.data <- genomic.data[genomic.data[,the.col] == data.chr,];
             the.chr <- paste(data.chroms[a.chr], "$", sep="");             
            cyto.data <- RCircos.Cyto[grep(the.chr, RCircos.Cyto$Chromosome),]

            if(max(in.data[,the.col+1]) > max(cyto.data[,3]) | 
                max(in.data[,the.col+2])> max(cyto.data[,3]))
            {
                message(paste(the.chr, max(in.data[,2]), max(in.data[,3])));
                stop("Location is outside of chromosome length.");  
            }
        }

        #   Make sure in genomic data all chromosome start positions 
        #   are smaller than their paired chromosome end positions
        #   =============================================================
        #
        for(a.row in 1:nrow(genomic.data))
        {
            if(genomic.data[a.row, the.col+1]>genomic.data[a.row, the.col+2]) 
            { 
                message("chromStart greater than chromEnd.\n"); 
                stop(paste("Row:", a.row, genomic.data[a.row, 2],  
                    genomic.data[a.row, 3]));
            }
        }
    }

    #   The validated data needs to be held for the RCircos session
    #   ===========================================================
    # 
    return (genomic.data);
}




#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
# 
#   18. Generate costomerized RCircos core components for multiple species
#       plot with all or part of chromosomes
#
#    Arguments:
#
#       cyto.info.list
#       List of multiple chromosome ideogram data
#
#       species
#       Character vector for prefix of chromosome names to identify species
#
#       chr.exclude
#       Chromosomes which should be excluded from dataset
#
#       tracks.inside and tracks.outside
#       Non-negative integer, number of data tracks inside and outside 
#       of chromosome ideogram
#
#    Return value:    None
#
#    Example:   cyto.info <- list(mouse.cyto, rat.cyto);
#               species.list <- c("m", "r");
#               RCircos.Multipl.Species.Core.Components(cyto.info, 
#                       species.list, NULL, 10, 0);
#
RCircos.Multiple.Species.Core.Components <- function(cyto.info.list, 
        species.list, chr.exclude=NULL, tracks.inside, tracks.outside)
{
    if(tracks.inside<0 || tracks.outside<0) 
    { stop("Track.inside and track.outside cannot be negative.\n"); }

    #   cyto.info.list and species must have same length
    #   ================================================
    #
    if(length(cyto.info.list) != length(species.list)) 
    { stop("Number of chromosome cytoband data and species must be same"); }

    #    Validte each chromosoem ideogram data then combine them as one
    #   ================================================
    #
    num.species <- length(cyto.info.list);
    for(a.cyto in 1:num.species)
    {
        cyto.info <- data.frame(cyto.info.list[a.cyto]);
        prefix <- species.list[a.cyto];

        cyto.info <- RCircos.Validate.Cyto.Info(cyto.info, NULL);
        cyto.info$Chromosome <- paste(prefix, cyto.info$Chromosome, sep="");

        if(a.cyto==1) { new.cyto.info <- cyto.info;
        } else { new.cyto.info <- rbind(new.cyto.info, cyto.info); }
    }

    #    Initialize RCircos core components
    #   =======================================================
    #
    RCircos.Initialize.Parameters(tracks.inside,tracks.outside);
    RCircos.Set.Cytoband.data(new.cyto.info);
    RCircos.Set.Base.Plot.Positions();
}


#   19. Combine and modify the chromosome names in multiple species datasets
#       to match the chromosomes in multiple species cytoband data
#   _______________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       data.list
#       List of genomic data from multiple species
#
#       species
#       Character vector for prefix of chromosome names to identify species
#
#       Note:   The order of each dataset in data list and species in 
#               species list must match.
#
#    Return value:  A data frame contain all data in the input data list
#                   with modified chromosome names 
#
#    Example:   data.list <- list(mouse.data, rat.data);
#               species.list <- c("m", "r");
#               dataset <- RCircos.Get.Multiple.Species.Dataset(data.list, 
#                                   species.list);
#
RCircos.Multiple.Species.Dataset <- function(data.list, species)
{
    #   Number of datasets and species must be same
    #   =======================================================
    #
    if(length(data.list) != length(species)) 
    { stop("Number of datasets and species must be same"); }

    #   Modify chromosome names in each dataset then combine 
    #   them as one
    #   =====================================================
    #
    num.data <- length(data.list);
    num.columns <- ncol(data.frame(data.list[1]));

    for(a.data in 1:num.data)
    {
        dataset <- data.frame(data.list[a.data]);
        prefix <- species[a.data];

        #    Number of columns of each dataset must be same
        #   ===============================================
        #
        if(ncol(dataset)!= num.columns) 
        { stop("Datasets have different columns."); }

        dataset[,1] <- paste(prefix, dataset[,1], sep="");
    
        if(a.data==1) { new.data <- dataset;
        } else { new.data <- rbind(new.data, dataset); }
    }

    return (new.data);
}


#   20. Ribbon link plot. Ribbons are wide link line between two 
#       chromosome segments.
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#    Arguments:
#
#       link.data
#       Data frame with paired genomic positions in each row
#
#       track.num
#       Number of the track from chromosome ideogram.
#
#   Return value:    None
#
#   Example:    RCircos.Ribbon.Plot(ribbon.data, 10, FALSE, FALSE)
#
RCircos.Ribbon.Plot <- function(ribbon.data, track.num, by.chromosome=FALSE, 
                            twist=FALSE)
{
    #   Arguments validate
    #   ==============================================
    #
    if(track.num<1) { stop("Track number cannot be smaller than 1.\n"); }
    if(by.chromosome != TRUE && by.chromosome != FALSE) 
    { stop("by.chromosome must be either TRUE or FALSE.\n"); }
    if(twist != TRUE && twist != FALSE) 
    { stop("twist must be either TRUE or FALSE.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Check chromosome names, Start, and End positions
    #   ================================================
    #
    ribbon.data <- RCircos.Validate.Genomic.Data(ribbon.data, plot.type="link");

    #   Plot position for link track.
    #   =============================================================
    #
    one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;
    track.out <- RCircos.Par$track.in.start - (track.num-1)*one.track;
    base.positions <- RCircos.Pos*track.out;

    #   Coordinates of the four conners of each ribbon (polygon) 
    #   ========================================================
    #
    data.points <- matrix(rep(0, nrow(ribbon.data)*4), ncol=4);
    for(a.link in 1:nrow(ribbon.data))
    {
        data.points[a.link, 1] <- RCircos.Data.Point(
                ribbon.data[a.link, 1], ribbon.data[a.link, 2]);
        data.points[a.link, 2] <- RCircos.Data.Point(
                ribbon.data[a.link, 1], ribbon.data[a.link, 3]);
        data.points[a.link, 3]<- RCircos.Data.Point(
                ribbon.data[a.link, 4], ribbon.data[a.link, 5]);
        data.points[a.link, 4]<- RCircos.Data.Point(
                ribbon.data[a.link, 4], ribbon.data[a.link, 6]);

        if(data.points[a.link, 1]==0 || data.points[a.link, 2]==0 ||
                data.points[a.link, 3]==0 || data.points[a.link, 4]==0)
        {  stop("Error in chromosome locations ...");  }
    }

    #    Ribbon colors
    #   ================================================
    #
    ribbon.colors <- RCircos.Get.Link.Colors(ribbon.data, by.chromosome);

    #   Draw each ribbon (polygon)
    #   ================================================
    #
    for(a.ribbon in 1:nrow(ribbon.data))
    {
        start.one <- data.points[a.ribbon, 1];
        end.one <- data.points[a.ribbon, 2];

        if(twist==FALSE) {
            start.two <- data.points[a.ribbon, 3];
            end.two <- data.points[a.ribbon, 4];
        } else {
            start.two <- data.points[a.ribbon, 4];
            end.two <- data.points[a.ribbon, 3];
        }

        P0 <- as.numeric(base.positions[end.one,]);
        P2 <- as.numeric(base.positions[start.two,]);
        line.one <- RCircos.Link.Line(P0, P2);

        P0 <- as.numeric(base.positions[end.two,]);
        P2 <- as.numeric(base.positions[start.one,]);
        line.two <- RCircos.Link.Line(P0, P2);

        polygon.x<- c(base.positions[start.one:end.one,1],
                        line.one$pos.x,
                        base.positions[start.two:end.two,1], 
                        line.two$pos.x );
        polygon.y<- c(base.positions[start.one:end.one,2], 
                        line.one$pos.y,
                        base.positions[start.two:end.two,2], 
                        line.two$pos.y );
        polygon(polygon.x, polygon.y, border=NA, 
                col=ribbon.colors[a.ribbon]);
    }
}


#   21. Erase one track or center area
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
# 
#   Arguments:
#
#       track.num
#       The ordinal numbers of track to be cleared (e.g., 2, 3, 4)
#
#       side
#       The position of track relative to chromosome ideogram, 
#       either "in" or "out"
#
#       to.center
#       Logical, TRUE for erase inner area including current 
#       track and FALSE for clear current track only
#
#    Returned values: None
#    Example:    RCircos.Clear.Track(track.num=2, side="in", to.center=FALSE);
#
RCircos.Clear.Track <- function(track.num, side, to.center=FALSE)
{
    if(track.num < 1) { stop("Track number cannot be smaller than 1.\n"); }
    if(to.center!=FALSE && to.center!=TRUE)
    {  stop("to.center must be either TRUE or FALSE.\n"); }

    side <- tolower(side);
    if(side != "in" && side!= "out") 
    { stop("side must be either in or out.\n"); }

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Plot position for current track. 
    #   ===================================================
    #
    locations <- RCircos.Track.Positions(side, track.num);

    if(side=="in") {
        out.pos <- locations[1];
        in.pos <- locations[2] - RCircos.Par$track.padding;
    } else {
        out.pos <- locations[1] + RCircos.Par$track.padding;
        in.pos <- locations[2];
    }

    #   Clear all inner area including current track 
    #   ============================================
    #
    if(side=="in" && to.center==TRUE)
    {
        polygon.x <- c(RCircos.Pos[,1]*out.pos, 
                        rep(0, nrow(RCircos.Pos)) );
        polygon.y <- c(RCircos.Pos[,2]*out.pos, 
                        rep(0, nrow(RCircos.Pos)) );        

    #    Clear current track area only
    #   ============================================
    #
    } else {
        polygon.x <- c(RCircos.Pos[,1]*out.pos, RCircos.Pos[,1]*in.pos);
        polygon.y <- c(RCircos.Pos[,2]*out.pos, RCircos.Pos[,2]*in.pos);
    }

    polygon(polygon.x, polygon.y, col="white", border="white");
}


#   22.    Reset chromosome ideogram plot data
#   ________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       chromIdeo
#       Data frame, an object of RCircos cytoband data returned from 
#       RCircos.Get.Plot.Ideogram() and modified for some items. 
#
#    Returned values: None
#
#   Example:    Reserved for advanced usage.
#               chromIdeo <- RCircos.Get.Plot.Ideogram();
#               chromIdeo$Location <- round(chromIdeo$Location*0.95);
#               RCircos.Reset.Plot.Ideogram(chromIdeo);
#
RCircos.Reset.Plot.Ideogram <- function(chromIdeo)
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Cytoband"]] <- chromIdeo;
}


#   23. Reset plot positions
#   _________________________________________________________________________
#    <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
#
#   Arguments:
#
#       plot.positions
#       A data frame with three columns for:
#           x and y coordinates for points which form a circle line 
#           degrees at each pioint the text should be rorated.
#
#   Returned values: None
#   Example:    Reserved for advanced usage.
#
RCircos.Reset.Plot.Positions <- function(plot.positions)
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Base.Position"]] <- plot.positions;
}


#   24. Draw ticks along chromosome highlight lines for a existing chromosome 
#       ideogram. Since in most cases the ticks are not needed, we add them
#       to existing image when necessary instead of drawing the chromosome 
#       ideogram with ticks at same time
#   ___________________________________________________________________________
#   <R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R><R>
##
#   Argument:    
#       tick.interval
#       distance between two ticks in MB base pairs
#
#   Return value:    None
#   Example:    RCircos.Tick.Plot(50);
#
RCircos.Add.Ideogram.Tick <- function(tick.interval)
{
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    #   Check if there is enough space for ticks and labels. From the  
    #   outside of chromosome ideogram, ticks start at highlight position 
    #   and take one track height, tick label takes three tracks, and  
    #   chromosome names use two tracks
    #   ==================================================================
    #
    ticks.height <- RCircos.Par$chr.ideog.pos + RCircos.Par$track.height*6;
    if(RCircos.Par$plot.radius < ticks.height ) 
    { 
        stop(paste0("There is no enough room to draw ticks.\n",
            "Please reset plot radius and redraw chromosome ideogram.\n")); 
    }

    #   Erase chromosome names
    #   ============================================================
    #
    chroms <- unique(RCircos.Cyto$Chromosome);
    polygon.out <- RCircos.Par$track.out.start;
    polygon.in  <- RCircos.Par$highlight.pos + 0.05;
    for(a.chr in 1:length(chroms))
    {
        the.chr  <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
        the.start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
        the.end   <- the.chr$Location[nrow(the.chr)];

        polygon.x<- c(RCircos.Pos[the.start:the.end,1]*polygon.out, 
                RCircos.Pos[the.end:the.start,1]*polygon.in);
        polygon.y<- c(RCircos.Pos[the.start:the.end,2]*polygon.out, 
                RCircos.Pos[the.end:the.start,2]*polygon.in);
        polygon(polygon.x, polygon.y, col="white", border="white");
    }

    #    Draw ticks and lables. Positions are calculated based on
    #    chromosome highlight positions
    #   ============================================================
    #
    start.pos <- RCircos.Par$highlight.pos;
    track.height <- RCircos.Par$track.height;

    in.pos   <- RCircos.Pos[, 1:2]*start.pos;
    out.pos  <- RCircos.Pos[, 1:2]*(start.pos + track.height);
    mid.pos  <- RCircos.Pos[, 1:2]*(start.pos + track.height/2);

    lab.pos  <- RCircos.Pos[, 1:2]*(start.pos + track.height*3);
    name.pos <- RCircos.Pos[, 1:2]*(start.pos + RCircos.Par$track.height*6);

    the.interval <- tick.interval*1000000;
    short.tick <- round(the.interval/RCircos.Par$base.per.unit, digits=0);
    long.tick  <- short.tick*2;

    chroms <- unique(RCircos.Cyto$Chromosome);
    for(a.chr in 1:length(chroms))
    {
        the.chr  <- RCircos.Cyto[RCircos.Cyto[,1]==chroms[a.chr],];
        chr.start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
        chr.end   <- the.chr$Location[nrow(the.chr)];

        #    Long ticks
        #   ==================================================
        #
        total.ticks <- ceiling((chr.end-chr.start)/long.tick);
        for(a.tick in 1:total.ticks)
        {
            #   Long ticks
            #   ================================================
            #
            tick.pos <- chr.start + (a.tick-1)*long.tick;
            if(tick.pos < chr.end)
            {
                lines(c(in.pos[tick.pos,1], out.pos[tick.pos,1]), 
                    c(in.pos[tick.pos,2], out.pos[tick.pos,2]),
                    col=the.chr$ChrColor[1]);

                #   Lables for long ticks
                #   ================================================
                # 
                lab.text <- paste0((a.tick-1)*tick.interval*2, "MB");
                text(lab.pos[tick.pos,1] , lab.pos[tick.pos,2], 
                        lab.text, cex=RCircos.Par$text.size,
                        srt=RCircos.Pos$degree[tick.pos]);
            }

            #    Short ticks
            #   =================================================
            #
            tick.pos <- tick.pos + short.tick;
            if(tick.pos < chr.end)
            {
                lines(c(in.pos[tick.pos,1], mid.pos[tick.pos,1]), 
                    c(in.pos[tick.pos,2], mid.pos[tick.pos,2]), 
                    col=the.chr$ChrColor[1]);
            }
        }

        #   Chromosome names
        #   =========================================================
        #
        name.loc  <- round((chr.end-chr.start)/2+chr.start, digits=0);
        text(name.pos[name.loc, 1], name.pos[name.loc, 2], 
                    sub("chr", "", chroms[a.chr]), 
                    srt=RCircos.Pos$degree[name.loc]);
    }
}



#   ======================================================================
#   <INTERNAL USE><INTERNAL USE><INTERNAL USE><INTERNAL USE><INTERNAL USE>


#   25. Initialize default global parameters for Circos plot. 
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
RCircos.Initialize.Parameters <- function(tracks.inside, tracks.outside)
{
    if(tracks.inside<0 ||tracks.outside<0) 
    { stop("Error: tracks.inside and tracks.outside cannot be negative.\n"); }

    radius.default <- 1.0;
    plot.radius.default <- 1.5;
    track.thick.default <- 0.12
    tracks.inside.default <- 4;

    #   Total of four data tracks are allowed inside of chromosome
    #   ideogram when radius.len is 1.0. Remains are for link lines. 
    #   If there are more than four tracks, increase the radius 
    #   length for extra plot area.
    #   ============================================================
    #
    more.inside <- tracks.inside - tracks.inside.default;
    if(more.inside > 0)
    { 
        more.area <- (more.inside + 1)* track.thick.default;
        radius.default <- radius.default + more.area;
        plot.radius.default <- plot.radius.default + more.area; 
    }

    #    If there will be data tracks outside of chromosome 
    #    ideogram, increase plot radius to get more room
    #   =======================================================
    #
    if(tracks.outside>0)
    { 
        more.area <- (tracks.outside + 1)* track.thick.default;
        plot.radius.default <- plot.radius.default + more.area; 
    }

    #    Set default plot parameters to a list
    #   ========================================================
    #
    plot.param <- list( 
            base.per.unit=30000,
            chrom.paddings=200,

            radius.len=radius.default,   
            chr.ideog.pos=radius.default+0.1,
            highlight.pos=radius.default+0.25,
            chr.name.pos=radius.default+0.4,

            plot.radius=plot.radius.default,
            track.out.start=radius.default+0.6, 
            track.in.start=radius.default+0.05, 

            chrom.width=0.1,
            track.padding=0.02, 
            track.height=0.1, 
 
            hist.width=1000,
            heatmap.width=1000,
            text.size=0.4,
            highlight.width=round(radius.default, digits=0),
            point.size=1,

            point.type=".", 
            text.color="black",
            heatmap.color="BlueWhiteRed",
            hist.color="red",
            line.color="black",
            scatter.color="black",
            tile.color="black",
            track.background="wheat",
            grid.line.color="gray",

            Bezier.point=1000,
            max.layers=5,
            sub.tracks=5);

    #    Put the plot parameter in RCircos environment
    #   =======================================================
    #
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.PlotPar"]] <- plot.param;
}


#   26. Set chromosome cytoband data for ideogram plot
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#    Arguments:    
#
#       cyto.band.info
#       A data frame which containing the chromosome ideogram data 
#       returned from function calls to RCircos.Validate.Cyto.Info()
#            
#    Returned values: None
#
RCircos.Set.Cytoband.data <- function(cyto.band.info)
{
    #   Reset colors for chromosome bands. Use yellow color for unknow 
    #   ==============================================================
    #
    stain2color <- as.character(cyto.band.info$Stain);
    band.color <- rep(colors()[652], length(stain2color));
    
    stains <- c("gneg", "acen", "stalk", "gvar", "gpos", "gpos100", 
                "gpos75", "gpos66", "gpos50", "gpos33", "gpos25");
    color.index <- c(1, 552, 615, 418, 24, 24, 193, 203, 213, 223, 233);

    for(a.stain in 1:length(stains))
    {
        bands <- which(stain2color==stains[a.stain]);
        if(length(bands)>0) 
        { band.color[bands] <- colors()[color.index[a.stain]]; }
    }
    cyto.band.info["BandColor"] <- band.color;

    #    Assign colors to chromosome highlight. There are total 50
    #    colors and the last 26 colors are reserved for future.
    #   ===========================================================
    #
    chrom.color <- c(552, 574, 645, 498, 450, 81, 26, 584, 524, 472,
            32, 57, 615, 635, 547, 254, 100, 72, 630, 589,
            8, 95, 568, 52);

    chrom2color <- as.character(cyto.band.info$Chromosome);
    chromosomes <- unique(chrom2color);

    #    In case of multiple ideogram plot, recycle the colors
    #   ===========================================================
    #
    num.chrom <- length(chromosomes);
    num.color <- length(chrom.color);
    if(num.chrom>num.color)
    {
        recycle.time <- floor(num.chrom/num.color);
        if(recycle.time>1) 
        { chrom.color <- rep(chrom.color, recycle.time); }

        remains <- num.chrom%%num.color
        if(remains > 0)
        {  chrom.color <- c(chrom.color, chrom.color[1:remains]); }
    }

    for(a.chr in 1:length(chromosomes))
    {
        rows <- which(chrom2color==chromosomes[a.chr]);
        if(length(rows)>0) 
        { chrom2color[rows] <- colors()[chrom.color[a.chr]]; }
    }
    cyto.band.info["ChrColor"] <- chrom2color;

    #    Total base pairs and relative length of each band
    #   ===========================================================
    #
    plot.par <- RCircos.Get.Plot.Parameters();

    cyto.band.info$ChromStart <- as.numeric(cyto.band.info$ChromStart);
    cyto.band.info$ChromEnd <- as.numeric(cyto.band.info$ChromEnd);

    band.len <- cyto.band.info$ChromEnd - cyto.band.info$ChromStart;
    cyto.band.info["Length"] <- band.len;
    cyto.band.info["Unit"]<- round(band.len/plot.par$base.per.unit, digits=0);

    #    Relative locations of each band in clockwise
    #   ===========================================================
    #
    Relative.Loc <- cyto.band.info$Unit;
    for(i in 2:length(Relative.Loc))
    { Relative.Loc[i] <- Relative.Loc[i] + Relative.Loc[i-1]; }
    cyto.band.info["Location"] <- Relative.Loc;

    if( plot.par$chrom.paddings>0) 
    {  
        chroms <- unique(cyto.band.info$Chromosome);
        chroms <- chroms[(chroms==chroms[1])==F];
        num.pad <-  plot.par$chrom.paddings;

        for(a.chr in 1:length(chroms))
        {
            index <- grep(paste(chroms[a.chr], "$", sep=""), 
                cyto.band.info$Chromosome);
            cyto.band.info$Location[index] <- num.pad + 
                cyto.band.info$Location[index];
            num.pad <- num.pad +  plot.par$chrom.paddings;
        }
    }

    #    Put the cyto.band.info data in RCircos environment
    #   ===========================================================
    #
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Cytoband"]] <- cyto.band.info;
}


#   27. Calculate the x and y coordinates for a circle line. These values 
#       will be used as base values to calculate the locations of plot
#       tracks, data points, and positions of chromosome band. Rotation 
#       degrees are also attached for text labeling at each point.
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#    Argument:    None
#    Return value:     None
#
RCircos.Set.Base.Plot.Positions <- function()
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #    Add one padding length to the last chromosome. Others
    #    are already included in RCircos.Cyto$Location
    #   ========================================================
    #
    total.points <- RCircos.Cyto$Location[nrow(RCircos.Cyto)] +  
                RCircos.Par$chrom.paddings;


    #    x and y coordinates for a circlar line with radius of 1,
    #    circumferance of 2*PI, and interval of 2PI/total.points
    #   =========================================================
    #
    interval <- 2*pi/total.points;
    base.val <- seq(0, 2*pi, interval);

    cor.x <- sin(base.val);
    cor.y <- cos(base.val);

    #    Degrees for text rotating at each posint
    #   =========================================================
    #
    degree <- rep(0, length(base.val));
    mid <- round((length(base.val)-1)/2, digits=0) + 1;

    total.points <- length(base.val);    
    degree[1:mid] <- 90 - (base.val[1:mid]*180/pi);
    degree[(mid+1):total.points] <- 270-(base.val[(mid+1):total.points]*180/pi);

    #    Put the plot postions data in RCircos environment
    #   ===========================================================
    #
    plot.postions <- data.frame(cor.x, cor.y, degree);

    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Base.Position"]] <- plot.postions;
}


#   29. Calculate plot coordinates for a chromosome position associated to a 
#       data point (e.g., the chromosome name and start position of a gene).
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#    Arguments:
#
#       chromosome   
#       Character vector, a chromosome name, e.g., "chr1"
#   
#       start
#       Non-negaive integer, start position of the gene on chromosomes
#
#   Return value:   An integer representing the index of x and y coordinates
#                   for a circular line
#  
RCircos.Data.Point <- function(chromosome, start)
{
    the.point <- 0;
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    #   Which band the start position is in
    #   ==========================================================
    #
    chrom.rows <- grep(paste("^", chromosome, "$", sep=""), 
            RCircos.Cyto$Chromosome);
    the.row <- which(RCircos.Cyto$ChromStart[chrom.rows] <= start &  
            RCircos.Cyto$ChromEnd[chrom.rows] >= start)[1];

    #    total length, units, and location of the band
    #   ==========================================================
    #
    band.length <- RCircos.Cyto$Length[chrom.rows[the.row]];
    band.units <- RCircos.Cyto$Unit[chrom.rows[the.row]];
    band.location <- RCircos.Cyto$Location[chrom.rows[the.row]];

    #    How far from the chromosome start the point is (by units)
    #   ==========================================================
    #
    the.bases <- start - RCircos.Cyto$ChromStart[chrom.rows[the.row]];
    the.units  <- the.bases/band.length*band.units;

    #    The exact point index of points for the circular line
    #   ==========================================================
    #
    the.point <- band.location - band.units + the.units;

    #    Return the index. Arguments are valudated outside of this 
    #    function so that there is no need to catch exception.
    #   ==========================================================
    #
    return (round(the.point, digits=0));
}


#   30. Convert genomic data to plot data by adding x and y coordinates for 
#       each row of a data set. A set of points for a circle is held in the 
#       RCircos session. We only need the index of the point for each
#       chromosome position
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#    Auguments:
#
#       genomic.data
#       A data frame contains genomic positions and data to be plotted
#
#        plot.type
#       Character vector, either "plot" or "link"                 
#
#    Returned value:    Updated genomic data (with a new column for Locations)
#
RCircos.Get.Plot.Data <- function(genomic.data, plot.type)
{
    #    Check chromosome names, chromStart, and chromEnd positions
    #   ==========================================================
    #
    genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type);

    #    Calculate the point index for each chromosome location
    #   ==========================================================
    #
    data.points <- rep(0, nrow(genomic.data));
    for(a.row in 1:nrow(genomic.data))
    {
        chromosome <- as.character(genomic.data[a.row, 1]);
        location <- round((genomic.data[a.row, 2] + 
                genomic.data[a.row, 3])/2, digits=0);
        data.points[a.row] <- RCircos.Data.Point(chromosome, location);
    }
    genomic.data["Location"] <- data.points;

    #    Sort the data by chromosome then start position
    #   ==========================================================
    #
    genomic.data <- genomic.data[order(genomic.data$Location),];

    #    The data needs to be held for the RCircos session
    #   ==========================================================
    #
    return (genomic.data);
}


#   31. Calculate inner and outer positions for a data track
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#  
#    Arguments:
#
#       side
#       Character vector, must be either "in" or "out".
#       
#       track.num
#       Non-negative integer, Number of the track from chromosome ideogram.
#
#    Return value:    the outer and iner positions of a track
#
RCircos.Track.Positions <- function(side, track.num)
{
    RCircos.Par <- RCircos.Get.Plot.Parameters()

    #    Track height.
    #   ==========================================================
    #
    one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;

    #    Positions based on side
    #   ==========================================================
    #
    side <- tolower(side);
    if(side=="in") 
    {
        out.pos <- RCircos.Par$track.in.start -(track.num-1)*one.track;
        in.pos  <- out.pos - RCircos.Par$track.height;
    } else if(side=="out"){
        in.pos <- RCircos.Par$track.out.start +(track.num-1)*one.track;
        out.pos <-in.pos + RCircos.Par$track.height;
    } else {  
        stop("Incorrect track location. It must be \"in\" or \"out\"."); 
    }

    if(out.pos <= 0 || in.pos <= 0) 
    { stop("Incorrect track location.\n"); }

    #    The position needs to be held for the RCircos session
    #   ======================================================
    #
    return (c(out.loc=out.pos, in.loc=in.pos));
}


#   32. Draw outline of one plot track. 
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#    Arguments:
#
#       out.pos
#       Non-negative numeric, scaling factor for outer position of the track
# 
#       in.pos
#       Non-negative numeric, scaling factor for inner position of the track 
#
#       num.layers
#       Non-negative integer, number of subtracks to plot, use 0 for no 
#       subtrack line
#
#       chrom.list
#       Character vector, of chromosome name list or NULL
#
#    Return value:    None
#
RCircos.Track.Outline <- function(out.pos, in.pos, num.layers, chrom.list)
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Subtrack height. Note: Some times there may have 
    #   more or less subtracks than defaul maximu layers
    #   =================================================
    #
    subtrack.height <- (out.pos-in.pos)/num.layers;
    chroms <- unique(RCircos.Cyto$Chromosome);

    #    Draw one or more but not all chromosomes
    #   =================================================
    #
    if(length(chrom.list)>0) { 
        if(sum(chrom.list %in% chroms)!= length(chrom.list))
        {  
            stop(paste("One or mroe chromosome is not",
                "in chromosome ideogram data.\n"));  
        }
        chroms <- chrom.list;
    } 

    for(a.chr in 1:length(chroms))
    {
        the.chr  <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
        start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
        end   <- the.chr$Location[nrow(the.chr)];

        polygon.x<- c(RCircos.Pos[start:end,1]*out.pos, 
                RCircos.Pos[end:start,1]*in.pos);
        polygon.y<- c(RCircos.Pos[start:end,2]*out.pos, 
                RCircos.Pos[end:start,2]*in.pos);
        polygon(polygon.x, polygon.y, col=RCircos.Par$track.background);

        if(num.layers>1) {
            for(a.line in 1:(num.layers-1))
            {
                height <- out.pos-a.line*subtrack.height;
                lines(RCircos.Pos[start:end,1]*height, 
                    RCircos.Pos[start:end,2]*height,
                    col=RCircos.Par$grid.line.color);
            }
        }
    }
}


#   33. In case of too many gene along one chromosome or a genomic region, 
#       the gene names may become very crowded so that we need reset plot 
#       positions to make the image more readable.
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#   Arguments:
#
#       genomic.data
#       A data frame with the four columns for chromosome name, start 
#       position, end position, and name of genes
#
#    Return value:  All or subset of input data frame with a new column 
#                   for plot positions.    
#
RCircos.Get.Gene.Label.Locations <- function(genomic.data)
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();

    units.default <- 1031904;
    width.default <- 5000;
    size.default  <- 0.4;

    #    Check chromosome names, Start, and End positions
    #   =================================================
    #
    genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, 
                            plot.type="plot");

    #    Calculate lable width (how many points it will be).
    #   =================================================
    #
    unit.factor <- units.default/sum(RCircos.Cyto$Unit);
    size.factor <- size.default/RCircos.Par$text.size;
    label.width <- (width.default/unit.factor)/size.factor;

    #    Get maximum number of lables for each chromosome. 
    #   =================================================
    #
    chromosomes <- unique(RCircos.Cyto$Chromosome);

    max.labels <- rep(0, length(chromosomes));
    start.loc <- rep(0, length(chromosomes));    
    end.loc  <- rep(0, length(chromosomes));    
    for(a.chr in 1:length(chromosomes))
    {
        index <- which(RCircos.Cyto$Chromosome==chromosomes[a.chr]);
        total.units <- sum(RCircos.Cyto$Unit[index]);
        max.labels[a.chr] <- floor(total.units/label.width);
        start.loc[a.chr] <- min(RCircos.Cyto$Location[index]);
        end.loc[a.chr] <- max(RCircos.Cyto$Location[index]);
    }

    #    Attach new column for label locations. 
    #   =================================================
    #
    label.data <- NULL;
    genomic.data["Label.Position"] <- genomic.data$Location;

    #    Reset label locations
    #   =================================================
    #
    message("\n");
    for(a.chr in 1:length(chromosomes))
    {
        index<- which(genomic.data[,1]==chromosomes[a.chr]);
        if(length(index)==0) { next; }

        #    If there too many gene labels, remove extra 
        #    genes for best visulization
        #   =================================================
        #
        if(length(index)>=max.labels[a.chr])
        {
            message(paste("Maximum lables for ", chromosomes[a.chr],
                "is", max.labels[a.chr],". " ));
            message("Extra ones are ignored.\n");

            the.chr <- genomic.data[index[1:max.labels[a.chr]],];
            the.chr <- the.chr[order(the.chr$Location),];

            for(a.gene in 1:nrow(the.chr))
            {
                new.loc  <- start.loc[a.chr] + (a.gene-1)* label.width;  
                the.chr$Label.Position[a.gene] <- new.loc;
            }
        } else {
            #    modify label locations if necessary
            #   =================================================
            #
            the.chr <- genomic.data[index,];
            the.chr <- the.chr[order(the.chr$Location),];
            all.gene <- nrow(the.chr);
            last.pos <- the.chr$Label.Position[1];
            for(a.gene in 1:all.gene)
            {
                curr.loc <-  the.chr$Label.Position[a.gene];

                len.needed <- curr.loc + 
                    (all.gene-a.gene)*label.width;
                if(len.needed>end.loc[a.chr]) 
                {  curr.loc  <- curr.loc - (len.needed-end.loc[a.chr]); }

                if(a.gene>1 && (curr.loc-last.pos) < label.width)
                { curr.loc  <- last.pos + label.width;  }
                
                the.chr$Label.Position[a.gene] <- curr.loc;
                last.pos <- curr.loc;                
            }
        }
        label.data <- rbind(label.data, the.chr); 
    }
    colnames(label.data) <- colnames(genomic.data);
    message("\n");

    #    The position needs to be held for the RCircos session
    #   =================================================
    #
    return (label.data);
}


#   34. Calculate x and y coordinated for a quandratic Bezier curve between 
#       two chromosome locations with the equation:  
#
#        B(t) = (1-t) ((1-t)P0 + tP1) + t((1-t)P1 + tP2)
#
#       where P0 is the start point, P2 is the end point, and P1 is the 
#       control point. Since we set P1 to (0,0), the equation become: 
#
#        B(t) =(1-t)^2*P0 + t^2*P2
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#   Arguments:
#
#       line.start
#       Numeric, x and y coordinates for the point where Bezier line starts
#
#       line.end
#       Numeric, x and y coordinates for the point where Bezier line ends
#
#    Return value:  a list contianing x and y coordinates for a quandratic 
#                   Bezier curve
#
RCircos.Link.Line <- function(line.start, line.end)
{    
    #   Set up the points for Bezure curve
    #   ==================================
    # 
    P0 <- line.start;
    P2 <- line.end;

    #   Calculate total number of points for the Bezuer curve
    #   =====================================================
    # 
    bc.point.num <- 1000;
    t <- seq(0, 1, 1/bc.point.num);

    #    Calculate the point values for Bezuer curve
    #   =====================================================
    # 
    link.x <- (1-t)^2*P0[1] + t^2*P2[1];
    link.y <- (1-t)^2*P0[2] + t^2*P2[2];

    #    Return the coordinates
    #   =====================================================
    # 
    return (list(pos.x=link.x, pos.y=link.y));
}


#   35. Create color map for heatmap plot
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#   Arguments:    
#
#       heatmap.color
#       Character vector, one of following:
#
#            BlueWhiteRed:      colors from blue to white then red
#            GreenWhiteRed:     colors from green to white then red
#            GreenYellowRed:    colors from green to yellow then red
#            GreenBlackRed:     colors from green to black then red
#            YellowToRed:       colors from yellow to red
#            BlackOnly:         default black colors
#
RCircos.Get.Heatmap.ColorScales <- function(heatmap.color)
{
    allOnes  <- seq(1, 1, length=256);
    allZeros <- seq(0, 0, length=256);
    one2zeor <- seq(1, 0, length=256);
    zero2one <- seq(0, 1, length=256);

    #   Blue, White, and Red
    #   ==========================================
    # 
    if(heatmap.color=="BlueWhiteRed") {

        RedRamp  <- rgb(allOnes, one2zeor, one2zeor);
        BlueRamp <- rgb(zero2one, zero2one, allOnes);    
        ColorRamp <- cbind(BlueRamp, RedRamp);

    #    Green, White, and Red
    #   ==========================================
    # 
    } else if (heatmap.color=="GreenWhiteRed") {

        RedRamp   <- rgb(allOnes, one2zeor, one2zeor);
        GreenRamp <- rgb(zero2one, allOnes, zero2one);
        ColorRamp <- cbind(GreenRamp, RedRamp);

    #   Green, Yellow, and Red
    #   ==========================================
    # 
    } else if (heatmap.color=="GreenYellowRed"){

        RedRamp <- rgb(allOnes, one2zeor, allZeros);
        GreenRamp <- rgb(zero2one, allOnes, allZeros);
        ColorRamp <- cbind(GreenRamp, RedRamp);

    #   Green, Black, and Red
    #   ==========================================
    # 
    } else if (heatmap.color=="GreenBlackRed"){

        RedRamp <- rgb( zero2one, allZeros, allZeros);
        GreenRamp <- rgb(allZeros, one2zeor, allZeros);
        ColorRamp <- cbind(GreenRamp, RedRamp);

    #    Yellow to Red
    #   ==========================================
    # 
    } else if (heatmap.color=="YellowToRed") {
        ColorRamp <- rgb(allOnes, one2zeor, allZeros);

    #    black only
    #   ==========================================
    # 
    } else { ColorRamp <- rgb(one2zeor, one2zeor, one2zeor); }

    return (ColorRamp);
}


#   36. Get plot color for each data point
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#   Arguments:    
#
#       plot.data
#       A data frame contining genomic data to be ploted
#       color
#       Character vector, default color name if it is not defined in plot data
#
#   Returned value: Character vector of color names or a R color
#
#
RCircos.Get.Plot.Colors <- function(plot.data, color)
{
        color.col <- grep("PlotColor", colnames(plot.data));

        if(length(color.col)==0)
        {  
                plot.colors <- rep(color, nrow(plot.data)); 
        } else if(length(color.col)==1) { 
                plot.colors <- as.character(plot.data[, color.col]); 
        } else { 
                stop("Incorrect plot colors defined in dataset.\n");
        }

        return (plot.colors);
}


#   37. Get plot color for each link line or ribbon
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
#
#   Arguments:    
#
#       plot.data
#       A data frame with paired genomic position data
#
#       by.chromosome
#       Logic, if TRUE, read color is used for links on same chromosome and 
#       blue color is used for linke between different chromosomes. If FALSE, 
#       user defined colors or rainbow colors will be used.
#
#   Returned value: Character vector of color names or vector of R colors
#
RCircos.Get.Link.Colors <- function(link.data, by.chromosome)
{
        red.color <- rgb(1,0,0, alpha=0.5);
        blue.color <- rgb(0, 0, 1, alpha=0.5);

        #   Default colors
        #   ============================================
        # 
        link.colors <- rep(blue.color, nrow(link.data));

        #   If by.chromosome is set to true, red color will be used 
        #   for links in same chromosome and blue color for links 
        #   between different chromosomes
        #   ============================================
        # 
        if(by.chromosome==TRUE) {

            for(a.row in 1:nrow(link.data))
            {
                if(link.data[a.row, 1]==link.data[a.row, 4]) 
                { link.colors[a.row] <- red.color; } 
            }

        #   If the plot color is provided in dataset, use it to
        #   replace the default one (rainbow)
        #   ============================================
        # 
        } else {
            color.col <- grep("PlotColor", colnames(link.data));

            if(length(color.col==1)) 
            {
                the.color <- as.character(link.data[, color.col]);
                for(a.row in 1:length(the.color))
                {
                    rgb.val <- as.vector(col2rgb(the.color[a.row]))/255;
                    link.colors[a.row] <- rgb(red=rgb.val[1], 
                        green=rgb.val[2], blue=rgb.val[3], alpha=0.5);
                }
            } else {
                for(a.row in 1:nrow(link.data))
                {
                    rgb.val <- as.vector(col2rgb(a.row+1))/255;
                    link.colors[a.row] <- rgb(red=rgb.val[1], 
                        green=rgb.val[2], blue=rgb.val[3], alpha=0.5);
                }
            }
        }

        return (link.colors);
}


#
#
#    Last revised on December 4, 2015
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>



