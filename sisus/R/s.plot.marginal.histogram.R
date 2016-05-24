s.plot.marginal.histogram <-
function# Plot marginal histograms
### internal function for sisus
(p.results
### internal variable
, p.results.name
### internal variable
, n.sources
### internal variable
, n.isotopes
### internal variable
, M.actual
### internal variable
, hist.bins
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
, SW
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  #library("splines")  # for density.pr()
  #library("locfdr")   # for density.pr()
  #library("fdrtool")  # for density.pr()

  color.histogram        = "steelblue3"; #gray50";
  color.density.smoother = "black";

  density.adjust=1; # for density smoother, multiplied by bandwidth

  for (i.plot in plot.format.list)
  {
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);

    #par(new = TRUE);
    #par(ask = TRUE);
    # Plotting histograms
    nbins = hist.bins;
    # for y-scaling and x-nbins
    par(mfrow=c(1, 1));
    max_y = 0;
    x.nbins = rep(0, n.sources);
    for (i in c(1:n.sources)){
      #x.nbins[i] = ceiling(nbins * (max(p.results[,i])-min(p.results[,i])));
      #r.max_y = hist(p.results[,i], x.nbins[i], plot=FALSE); max_y = max(r.max_y$density,max_y); # calculate max_y for equal height of histograms

        r.max_y = hist(p.results[,i], nbins, plot=FALSE); max_y = max(r.max_y$density,max_y); # calculate max_y for equal height of histograms
    }

    #par(oma=c(5,4,4,4), mfrow=c(n.sources,1), bg="white")
    #par(mfrow=c(n.sources,1));
    par(mfrow=c(n.sources,1), mar=c(.5,5,.5,2), oma=c(7,4,5,4));  # mar allows the histograms to touch top-bottom c(bot,lef,top,rig)
    for (i in 1:n.sources) {
      #omar =
      #par(mar=c(2,4,0,0));
      #hist(p.results[,i], nbins, freq = FALSE, main=NULL, xlab = NULL, ylab=paste(names.sources[i]), xlim = c(0,1), ylim = c(0,max_y), axes = TRUE, plot = TRUE, col = "black", labels = FALSE);
      #p.results.density = density.pr(p.results[,i], plot = TRUE, ncells = x.nbins[i], main=NULL, xlab = NULL, ylab=paste(names.sources[i]), xlim = c(0,1), ylim = c(0,max_y), axes = TRUE, col = "black", labels = FALSE);

      # p.results.density = density.pr(p.results[,i], plot = TRUE, ncells = nbins, main=NULL, xlab = NULL, ylab=paste(names.sources[i]), xlim = c(0,1), ylim = c(0,max_y), axes = TRUE, col = "black", labels = FALSE, lwd=2);

 # 4/29/2007 10:40AM removed since density.pr() is no longer available
 ##     p.results.density = density.pr(p.results[,i], plot = TRUE, ncells = nbins, main=NULL, xlab = NULL, ylab=paste(names.sources[i]), xlim = c(0,1), ylim = c(0,max_y), axes = TRUE, col = "black", labels = FALSE, lwd=2, xpd=NA );
      # density.pr() no longer works, use density() to smooth histogram
      r.max_y = hist(p.results[,i], nbins, freq=FALSE, plot=TRUE, main=NULL, xlab = NULL, ylab=paste(names.sources[i]), xlim = c(0,1), xaxp = c(0,1,10), ylim = c(0,max_y), axes = TRUE, col = color.histogram, labels = FALSE);
      if (SW$hist.density.sw == 1) {
        lines(density(p.results[,i], kernel="optcosine", from=max(0,min(p.results[,i])), to=min(1,max(p.results[,i])), adjust=density.adjust), lwd=2, xpd=NA, col = color.density.smoother ); # density smooth histogram
      }
          ## possible kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")

      #par(xpd=TRUE);
      #par(omar);
    }
    # box("outer", col="black")  # draw a box around the entire figure
    mtext(paste("Marginal Histograms of feasible source contribution proportions"), side=1, line=3, outer=TRUE);  # bottom
    mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=5, outer=TRUE);  # bottom (line specifies)
    #mtext(paste("Source contribution (%)"), side=1, line=1, outer=TRUE);                         # bottom
    #mtext(paste("Stable Isotope Sourcing using Sampling (SISUS)"), side=1, line=3, outer=TRUE);  # bottom (line specifies)
    mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                             # left
    mtext(paste(p.results.name,": ", paste(names.mixtures.indy)), cex = 1.2, side=3, line=2, outer=TRUE);          # top
    mtext(paste(analysis.name), cex = 0.8, side=3, line=.5, outer=TRUE);                                      # top
    mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=0, outer=TRUE);  # right
    #par(ask = FALSE);

    #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
  } # i.plot plotting loop

  ### internal variable
}
