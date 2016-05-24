mcmc.diagnostics <-
function# MCMC Diagnostics
### internal function for sisus
(p.biomass.sam
### internal variable
, p.results.name
### internal variable
, names.mixtures.indy
### internal variable
, n.sources
### internal variable
, names.sources
### internal variable
, n.isotopes
### internal variable
, names.isotopes
### internal variable
, M.actual
### internal variable
, output.mcmc.diagnostics.filename
### internal variable
, filename.prefix
### internal variable
, analysis.name
### internal variable
, plot.format.list
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

#### cross corr plots are not plotting correctly
#### Geweke plot has cutoffs extending beyond plot box

#     Provide a collection of diagnostic output with interpretations and suggestions if failed
#       R coda package: (http://cran.r-project.org/src/contrib/Descriptions/coda.html)

  # attached via DESCRIPTION # library("coda");
  # attached via DESCRIPTION # library("stats");  # for capture.output()

  alpha.sig = 0.05;
  acf.max = 0.05;

  # as.mcmc() to create as an mcmc object to use by the functions below
  #mcmc.p = as.mcmc(p.biomass.sam);
  mcmc.p = mcmc(data=p.biomass.sam, start=1, end=M.actual, thin=1);

  output.mcmc =
    rbind(
      paste("SISUS: MCMC Diagnostics of ", paste(sprintf("%s for %s", p.results.name, names.mixtures.indy)))
     ,paste("       ", paste(sprintf("%s", analysis.name)), sep="")
    );
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);


  ########################################
  # make source names a common length, either by shortening or by padding with spaces
  var.len = 8; # column length for output
  nchar.names.sources = nchar(names.sources);
  names.sources.padded = names.sources;
  for (i in seq(1,n.sources)) {
    if (nchar.names.sources[i] < var.len) {
      names.sources.padded[i] = paste(paste(rep(" ", var.len-nchar.names.sources[i]-1),collapse=""), names.sources[i]);
    }
  }

  ################################################################################
  ########################################
  # autocorr() for autocorrelation
  plot.filename = paste("plot_mcmc_diagnostics_", names.mixtures.indy, "_autocorrelations", sep = ""); plot.filename = filename.clean(plot.filename);
  max.lag = 25;
  output.mcmc =
    rbind(
      paste(""),paste("")
     ,paste("================================================================================")
     ,paste("Autocorrelations from lag 1 to lag ", paste(sprintf("%s", max.lag)), sep="")
     ,paste("  These should all be close to 0, if not then increase skip parameter.")
     ,paste("  See plot: ", plot.filename)
     ,paste("")
     #,paste("Lag ", paste(sprintf("%s", substr(names.sources.padded, 1, var.len)), collapse="   "))
    );
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);
  r.autocorr = round(autocorr(mcmc.p[,], lags=seq(1,max.lag), relative=TRUE)[,,1], digits=5);
  #output.mcmc = r.autocorr;
  #output.mcmc = round(r.autocorr, digits=5);

  #for (i.lag in seq(1,max.lag)) {
  #write(paste(sprintf("%3d ",i.lag), paste(sprintf("% 1.5f",output.mcmc[i.lag,])     , collapse="   "))
  #      , file = output.mcmc.diagnostics.filename, ncolumns = n.sources, append = TRUE, sep = "   ");
  #
  #}

  capture.output(r.autocorr, file = output.mcmc.diagnostics.filename, append = TRUE);

  status.acf = NULL;
  if (sum(abs(r.autocorr[1,]) < acf.max) == n.sources) {  # if they all pass
    output.mcmc = rbind(paste(""),paste("  * pass *"),paste(""));
    status.acf = 0;
  } else {
    output.mcmc = rbind(paste(""),paste("  !!! FAIL !!!  Increase skip"),paste(""));
    status.acf = 1;
  }
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);



  ########################################
  # autocorr.plot() to plot matrix of autocorrelations
  for (i.plot in plot.format.list)
  {
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);

    par(mfrow=c(n.sources,1), mar=c(2,5,2,2), oma=c(7,4,5,4));  # mar allows the histograms to touch top-bottom c(bot,lef,top,rig)

    ## Create plot
    autocorr.plot(mcmc.p, lag.max=max.lag, auto.layout=FALSE);

    # box("outer", col="black")  # draw a box around the entire figure
    mtext(paste("MCMC Diagnostics: Autocorrelations"), side=1, line=3, outer=TRUE);  # bottom
    mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=5, outer=TRUE);  # bottom (line specifies)
    mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                             # left
    mtext(paste(p.results.name,": ", paste(names.mixtures.indy)), cex = 1.2, side=3, line=2, outer=TRUE);          # top
    mtext(paste(analysis.name), cex = 0.8, side=3, line=.5, outer=TRUE);                                      # top
    mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=0, outer=TRUE);  # right
    #par(ask = FALSE);

    #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
  } # i.plot plotting loop

# Removed crosscorr since it may not provide useful information for diagnostics because of the system constraints
# Yet, it may be useful as supplimental information to the scatterplot matrix
#  ################################################################################
#  ########################################
#  # crosscorr() for correlations between variables
#  plot.filename = paste("plot_mcmc_diagnostics_", names.mixtures.indy, "_crosscorrelations", sep = ""); plot.filename = filename.clean(plot.filename);
#  output.mcmc =
#    rbind(
#      paste(""),paste("")
#     ,paste("================================================================================")
#     ,paste("Cross Correlations (not diagnostic, information only)")
#     ,paste("  Indicates how source solutions correlate (reference Biomass scatterplot).")
#     ,paste("  See plot: ", plot.filename)
#     ,paste("")
#     #,paste("          ", paste(sprintf("%s", substr(names.sources.padded, 1, var.len)), collapse="   "))
#    );
#  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);
#  r.crosscorr = crosscorr(mcmc.p);
#  #output.mcmc = r.crosscorr;
#  #output.mcmc = round(r.crosscorr, digits=5);
#  #write(format(output.mcmc, digits=5), file = output.mcmc.diagnostics.filename, ncolumns = n.sources, append = TRUE, sep = "   ");
#
#  #for (i.source in seq(1,n.sources)) {
#  #write(paste(sprintf("%s  ",names.sources.padded[i.source]), paste(sprintf("% 1.5f",output.mcmc[i.source,])     , collapse="   "))
#  #      , file = output.mcmc.diagnostics.filename, ncolumns = n.sources, append = TRUE, sep = "   ");
#  #
#  #}
#
#  capture.output(r.crosscorr, file = output.mcmc.diagnostics.filename, append = TRUE);


  ## ########################################
  ## # levelplot(line[[2]]) also for crosscorrelation plot
  ## for (i.plot in plot.format.list)
  ## {
  ##   s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);
  ##
  ##   #par(mar=c(2,5,2,2), oma=c(7,4,5,4));  # mar allows the histograms to touch top-bottom c(bot,lef,top,rig)
  ##
  ##   ## Create plot
  ##   #levelplot(mcmc.p);
  ##   levelplot(mcmc.p, main=paste(p.results.name,": ", paste(names.mixtures.indy)), sub=paste("MCMC Diagnostics: Cross Correlations"));
  ##
  ##   # box("outer", col="black")  # draw a box around the entire figure
  ##  # mtext(paste("MCMC Diagnostics: Cross Correlations"), side=1, line=3, outer=TRUE);  # bottom
  ##  # mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=5, outer=TRUE);  # bottom (line specifies)
  ##  # mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                             # left
  ##  # mtext(paste(p.results.name,": ", paste(names.mixtures.indy)), cex = 1.2, side=3, line=2, outer=TRUE);          # top
  ##  # mtext(paste(analysis.name), cex = 0.8, side=3, line=.5, outer=TRUE);                                      # top
  ##  # mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=0, outer=TRUE);  # right
  ##   #par(ask = FALSE);
  ##
  ##   #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
  ##   s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
  ## } # i.plot plotting loop


  ################################################################################
  ########################################
  # geweke.diag() for testing mean of first 0.1 same as last 0.5
  plot.filename = paste("plot_mcmc_diagnostics_", names.mixtures.indy, "_geweke", sep = ""); plot.filename = filename.clean(plot.filename);
  geweke.frac1=0.1; geweke.frac2=0.5;
  output.mcmc =
    rbind(
      paste(""),paste("")
     ,paste("================================================================================")
     ,paste("Geweke's convergence diagnostic")
     ,paste("  Testing whether mean of first ", geweke.frac1, " of samples is different from last ", geweke.frac2)
     ,paste("    Indicates whether the first and last part of a sample from the Markov chain are drawn from the same distribution (z-test)")
     ,paste("    if reject, then see plot", plot.filename, "to suggest length of burnin = [plot value in cutoffs]*skip")
     ,paste("")
     ,paste("        ", paste(sprintf("%s", substr(names.sources.padded, 1, var.len)), collapse="   "))
    );
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);

  r.geweke.diag = geweke.diag(mcmc.p, frac1=geweke.frac1, frac2=geweke.frac2);
  r.geweke.diag.pvalue = 2*pnorm(abs(r.geweke.diag$z), mean=0, sd=1, lower.tail = FALSE);
  output.mcmc =
    rbind(
      paste("z-score ", paste(sprintf("% 1.5f",r.geweke.diag$z     )     , collapse="   "))
     ,paste("p-value ", paste(sprintf("% 1.5f",r.geweke.diag.pvalue)     , collapse="   "))
    )
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);


  status.geweke = NULL;
  if (sum(r.geweke.diag.pvalue >= alpha.sig) == n.sources) {  # if they all pass
    output.mcmc = rbind(paste(""),paste("  * pass *"),paste(""));
    status.geweke = 0;
  } else {
    output.mcmc = rbind(paste(""),paste("  !!! FAIL !!!  Increase burnin"),paste(""));
    status.geweke = 1;
  }
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);

  #capture.output(r.geweke.diag, file = output.mcmc.diagnostics.filename, append = TRUE);


  ########################################
  for (i.plot in plot.format.list)
  {
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);

    par(mfrow=c(n.sources,1), mar=c(2,5,2,2), oma=c(7,4,5,4));  # mar allows the histograms to touch top-bottom c(bot,lef,top,rig)

    ## Create plot
    geweke.plot(mcmc.p, frac1=geweke.frac1, frac2=geweke.frac2, nbins=20, pvalue=alpha.sig, auto.layout=FALSE);

    # box("outer", col="black")  # draw a box around the entire figure
    mtext(paste("MCMC Diagnostics: Geweke's convergence diagnostic"), side=1, line=3, outer=TRUE);  # bottom
    mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=5, outer=TRUE);  # bottom (line specifies)
    mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                             # left
    mtext(paste(p.results.name,": ", paste(names.mixtures.indy)), cex = 1.2, side=3, line=2, outer=TRUE);          # top
    mtext(paste(analysis.name), cex = 0.8, side=3, line=.5, outer=TRUE);                                      # top
    mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=0, outer=TRUE);  # right
    #par(ask = FALSE);

    #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
  } # i.plot plotting loop


  ################################################################################
  ########################################
  # heidel.diag() for relative accuracy for estimating the mean
  heidel.eps = 0.1;
  output.mcmc =
    rbind(
      paste(""),paste("")
     ,paste("================================================================================")
     ,paste("Heidelberger and Welch's convergence diagnostic")
     ,paste("  Convergence test uses the Cramer-von-Mises statistic to test the")
     ,paste("    null hypothesis that the sampled values come from a stationary distribution.")
     ,paste("  Halfwidth test is a run length control diagnostic based on a ")
     ,paste("    criterion of relative accuracy for the estimate of the mean.")
     ,paste("    Test corresponds to a relative accuracy of two significant digits.")
     ,paste("")
     #,paste("        ", paste(sprintf("%s", substr(names.sources.padded, 1, var.len)), collapse="   "))
    );
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);

  r.heidel.diag = heidel.diag(mcmc.p, eps=heidel.eps, pvalue=alpha.sig);

  capture.output(r.heidel.diag, file = output.mcmc.diagnostics.filename, append = TRUE);


  r.heidel.diag.conv = r.heidel.diag[seq(1,n.sources)];
  r.heidel.diag.half = r.heidel.diag[n.sources*3+seq(1,n.sources)];
  status.heidel = NULL;
  if (sum(r.heidel.diag.conv + r.heidel.diag.half) == 2*n.sources) {  # if they all pass
    output.mcmc = rbind(paste(""),paste("  * pass *"),paste(""));
    status.heidel = 0;
  } else {
    output.mcmc = rbind(paste(""),paste("  !!! FAIL !!!  Increase burnin"),paste(""));
    status.heidel = 1;
  }
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);


  ################################################################################
  ########################################
  # plot() of mcmc object does trace and density plots for variables
  #plot.filename = paste("plot_mcmc_diagnostics_", names.mixtures.indy, "_type", sep = ""); plot.filename = filename.clean(plot.filename);
  #plot(mcmc.p);

  ################################################################################
  ########################################
  # raftery.diag() determines sample size required for estimating quantiles with specified precision
  output.mcmc =
    rbind(
      paste(""),paste("")
     ,paste("================================================================================")
     ,paste("Raftery and Lewis's diagnostic")
     ,paste("  Run length control diagnostic based on a criterion of accuracy of estimation of the quantile q.")
     ,paste("  The number of iterations required to estimate ")
     ,paste("    the quantile q to within an accuracy of +/- r with probability p is calculated.")
     ,paste("  The number of `burnin' iterations to be discarded at the beginning of the chain is also calculated.")
     ,paste("")
     #,paste("        ", paste(sprintf("%s", substr(names.sources.padded, 1, var.len)), collapse="   "))
    );
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);

  r.raftery.diag = raftery.diag(mcmc.p, q=0.01, r=0.005, s=1-alpha.sig, converge.eps=0.001);

  capture.output(r.raftery.diag, file = output.mcmc.diagnostics.filename, append = TRUE);

  status.raftery = NULL;
  if (!(r.raftery.diag$resmatrix[1]=="Error")) {
    if (sum(r.raftery.diag$resmatrix[,3] <= M.actual) == n.sources) {  # if they all pass
      output.mcmc = rbind(paste("  * pass *"),paste(""));
      status.raftery = 0;
    } else {
      output.mcmc = rbind(paste("  !!! FAIL !!!  Increase sample size M, at least the number suggested by Nmin"),paste(""));
      status.raftery = 1;
    };
  } else {
    output.mcmc = rbind(paste("  !!! FAIL !!!  Increase sample size M, at least the number suggested by Nmin = ", r.raftery.diag$resmatrix[2]),paste(""));
    status.raftery = 1;
  };
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);


  ################################################################################
  ########################################
  # traceplot() for chains
  plot.filename = paste("plot_mcmc_diagnostics_", names.mixtures.indy, "_traceplot", sep = ""); plot.filename = filename.clean(plot.filename);
  output.mcmc =
    rbind(
      paste(""),paste("")
     ,paste("================================================================================")
     ,paste("Trace plot of MCMC output")
     ,paste("  Displays a plot of iterations vs. sampled values for each source in the chain.")
     ,paste("  Should be evenly scattered throughout for each source")
     ,paste("  See plot: ", plot.filename)
     ,paste("")
    );
  write(output.mcmc, file = output.mcmc.diagnostics.filename, append = TRUE);

  ########################################
  for (i.plot in plot.format.list)
  {
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);

    par(mfrow=c(n.sources,1), mar=c(2,5,2,2), oma=c(7,4,5,4));  # mar allows the histograms to touch top-bottom c(bot,lef,top,rig)

    ## Create plot
    traceplot(mcmc.p);

    # box("outer", col="black")  # draw a box around the entire figure
    mtext(paste("MCMC Diagnostics: Trace plots of chain"), side=1, line=3, outer=TRUE);  # bottom
    mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=5, outer=TRUE);  # bottom (line specifies)
    mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                             # left
    mtext(paste(p.results.name,": ", paste(names.mixtures.indy)), cex = 1.2, side=3, line=2, outer=TRUE);          # top
    mtext(paste(analysis.name), cex = 0.8, side=3, line=.5, outer=TRUE);                                      # top
    mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=0, outer=TRUE);  # right
    #par(ask = FALSE);

    #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
    s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
  } # i.plot plotting loop



  ################################################################################
  ########################################
  # HPDinterval() for giving most likely credible intervals
  #HPDinterval(mcmc.p, prob=0.95)


  mcmc.diag.status = new.env();  # create a list to return status of each test

  mcmc.diag.status$acf     = status.acf    ;
  mcmc.diag.status$geweke  = status.geweke ;
  mcmc.diag.status$heidel  = status.heidel ;
  mcmc.diag.status$raftery = status.raftery;

  return( as.list(mcmc.diag.status) );


  ### internal variable
}
