# This function originates from the GLAD package (GPL).  We should
# rewrite it from scratch to avoid license issues. /HB 2010-02-19
setMethodS3("drawCytoband2", "default", function(cytoband, chromosome=1, y=-1, labels=TRUE, height=1, colCytoBand=c("white", "darkblue"), colCentro="red", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  seqPalette <- NULL;
  # Try to use the palette defined by GLAD
  if (isPackageInstalled("GLAD")) {
    # Note that although GLAD is installed we may still get:
    # > GLAD::myPalette
    # Error in library.dynam(lib, package, package.lib) :
    # DLL 'GLAD' not found: maybe not installed for this architecture?

    # In order to minimize the impact of that, we will indeed to load
    # GLAD although we're only interested in importing GLAD::myPalette().
    # We could unload GLAD when exiting this function (iff it was loaded
    # by us), but since drawCytoband2() is called for each sample and
    # chromosome, that would generate a huge number of package loads.
    # /HB 2010-12-07

    tryCatch({
      requireWithMemory("GLAD") || throw("Package not loaded: GLAD");
      # WORKAROUND: Since GLAD is not using packageStartupMessage()
      # but cat() in .onLoad(), there will be a long message printed
      # even when using GLAD::<fcn>. /HB 2010-02-19
      dummy <- capture.output({
        fcn <- GLAD::myPalette;
      });
      seqPalette <- function(from, to, length.out, ...) {
        fcn(low=from, high=to, k=length.out, ...);
      } # seqPalette()
    }, error=function(ex) {});
  }

  # Fallback...
  if (is.null(seqPalette)) {
    seqPalette <- function(from, to, length.out, ...) {
      # Argument 'from':
      if (is.character(from)) {
        from <- col2rgb(from) / 255;
      }
      # Argument 'to':
      if (is.character(to)) {
        to <- col2rgb(to) / 255;
      }

      # Generate sequence of RGB vectors
      rgbData <- mapply(from, to, FUN=function(from, to) {
        seq(from=from, to=to, length.out=length.out);
      });

      # Translate to colors
      rgb(rgbData);
    } # seqPalette()
  }

  opar <- par(xpd=NA);
  on.exit(par(opar));


  # Nothing todo?
  if (nrow(cytoband) == 0) {
    return();
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cytoband colors
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  color <- unique(cytoband$Color);
  pal <- seqPalette(from=colCytoBand[1], to=colCytoBand[2],
                                                 length.out=length(color));

  info <- data.frame(Color=color, ColorName=I(pal));
  cytoband <- merge(cytoband, info, by="Color");
  # Not needed anymore
  info <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract cytoband information for current chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keep <- which(cytoband$Chromosome == chromosome);
  cytoband <- cytoband[keep, ,drop=FALSE];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cytoband positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  CytoPos <- 0.5 * (cytoband$Start + cytoband$End);
  CytoLength <- (cytoband$End - cytoband$Start);
  NbCyto <- length(cytoband[, 1]);
  HeightPlot <- rep(height, NbCyto);
  sizeCyto <- matrix(c(CytoLength, HeightPlot), nrow=NbCyto, ncol=2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Draw cytobands
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  y0 <- min(unique(y));
  yC <- y0+height/2;
  y1 <- y0+height;
  symbols(x=CytoPos, y=rep(yC, NbCyto), rectangles=sizeCyto,
      inches=FALSE, bg=cytoband$ColorName, add=TRUE, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Highlight the centromere
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The inverted arrow indicating where the centromere is.
  idxs <- which(cytoband$Centro == 1);
  centroPos <- min(cytoband$End[idxs]);
  arrows(centroPos, y0, centroPos, y1, col=colCentro, code=2, angle=120,
                                                               length=0.1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Labels, e.g. 20q12
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (labels) {
    labels <- paste(cytoband$Chromosome, cytoband$Band, sep="");
#    axis(side=3, at=CytoPos, labels=labels, las=2);
    dy <- par("cxy")[2];
    text(x=CytoPos, y=y1+dy/2, labels=labels, srt=90, adj=c(0,0.5));
  }
}, private=TRUE)



############################################################################
# HISTORY:
# 2010-12-07
# o Now drawCytoband2() utilizes requireWithMemory() to minimize the
#   annoyances when GLAD fails to load.
# 2010-12-02
# o BUG FIX: drawCytoband2() would throw an error if argument 'cytoband'
#   was an empty data frame.  Now it returns quietly.
# 2010-10-13
# o ROBUSTNESS/BUG FIX: The internal drawCytoband2() used to annotate
#   chromosomal plots with cytobands tries to utilize GLAD package,
#   if available.  However, even when GLAD is installed it may still be
#   broken to missing dynamic libraries, e.g. 'Error in library.dynam(lib,
#   package, package.lib) : DLL 'GLAD' not found: maybe not installed for
#   this architecture?'.  We now avoid this too.
# 2010-02-19
# o Now suppressing the .onLoad() message of GLAD.
# o Extracted from GLAD.EXT.R.
# o Now the default drawCytoband2() can work without GLAD, if not installed.
# 2009-05-10
# o Moved to aroma.core v1.0.6.  Source files: profileCGH.drawCnRegions.R
#   and profileCGH.drawCytoband.R.
# 2007-08-22
# o Update plotProfile2() to utilizes drawCnRegions().
# 2007-06-11
# o Added explicit call to GLAD::myPalette() to please R CMD check R v2.6.0.
# 2007-01-03
# o Made the highlighting "arrow" for the centromere smaller.
# 2006-12-20
# o It is now possible to specify 'xlim' as well as 'ylim'.
# o Reimplemented, because the cytoband was not displayed correctly.
############################################################################
