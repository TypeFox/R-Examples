s.plot.scatterplot.sample <-
function# Scatterplot matrix
### internal function for sisus
(p.results
### internal variable
, p.results.name
### internal variable
, n.sources
### internal variable
, n.isotopes
### internal variable
, names.isotopes
### internal variable
, names.mixtures.indy
### internal variable
, names.sources
### internal variable
, analysis.name
### internal variable
, filename.prefix
### internal variable
, plot.filename
### internal variable
, plot.format.list
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # Scatterplots with smoothed densities color representation (http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=139)
  # library("geneplotter"); # will load "annotate", "Biobase", "tools"  REMOVED  3/5/2014
  # attached via DESCRIPTION # library("MASS")  # for kde2d() density plot
  # attached via DESCRIPTION # library("RColorBrewer");

  for (i.plot in plot.format.list)
  {
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);

    axis.limits = c(min(p.results), max(p.results));
    #par(new = TRUE);
    #par(bg="white")
    #par(oma=c(5,4,4,4)*2, bg="white")
    #par(mar=c(2,4,0,0)*2);
    #par(ask = TRUE)
    pairs(p.results,
          labels = names.sources,
          oma = c(5,4,4,5)*2, mar = c(4,4,2,2)*2,
          #upper.panel = function(...) { par(new=TRUE); smoothScatter(..., nrpoints=Inf, xaxt = "n", yaxt = "n"); },
          upper.panel = function(x,y, ...) { par(new=TRUE); plot(x, y, type="p", pch=19, cex=.001, xaxt = "n", yaxt = "n"); },
          lower.panel = function(x,y, ...) { par(new=TRUE); image(kde2d(x, y, h=0.025, n=100, lims=c(0,1,0,1)), col = colorRampPalette(brewer.pal(9,"Blues"))(100), add=TRUE) },
                                                     # smoothScatter(..., nrpoints=0  , bandwidth=c(0.01,0.01), xaxt = "n", yaxt = "n"); },
          #upper.panel = function(...) { par(new=TRUE); smoothScatter(..., nrpoints=0); },
          #lower.panel = function(...) { par(new=TRUE); colors = densCols(...); plot(..., col=colors, pch=200); }
          #lower.panel = function(...) { par(new=TRUE); plot(..., type="p", pch=19, cex=.001); },
          #lower.panel = function(...) { par(new=TRUE); hist2d(..., nbins=20, col=c("white", heat.colors(16)) ); }
          diag.panel = diag.panel.hist,
          xlim = axis.limits,
          ylim = axis.limits,
          gap = 0
          );
    #par(xpd=TRUE);

    # box("outer", col="black")  # draw a box around the entire figure
    mtext(paste("Scatterplot Matrix"), side=1, line=1, outer=TRUE);                         # bottom
    mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=3, outer=TRUE);  # bottom (line specifies)
    mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                             # left
    mtext(paste(p.results.name,": ", paste(names.mixtures.indy)), cex = 1.2, side=3, line=1, outer=TRUE);          # top
    #mtext(paste(paste(names.mixtures.indy, ": ", tit)), cex = 1.2, side=3, line=1, outer=TRUE);          # top
    mtext(paste(analysis.name), cex = 0.8, side=3, line=0, outer=TRUE);                         # top   # line=-.5 in other plots
    mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=0, outer=TRUE);  # right
    #par(ask = FALSE)

    #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
  } # plotting loop

  ### internal variable
}
