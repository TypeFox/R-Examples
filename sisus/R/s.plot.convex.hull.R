s.plot.convex.hull <-
function# Convex Hull plot
### internal function for sisus
(ch.isotopes.sources
### internal variable
, ch.isotopes.mixtures
### internal variable
, ch.n.mixtures
### internal variable
, n.sources
### internal variable
, n.isotopes
### internal variable
, names.sources
### internal variable
, names.isotopes
### internal variable
, ch.names.mixtures
### internal variable
, title.mixture
### internal variable
, ch.isotope.mvn.sample
### internal variable
, ch.concentration.sources
### internal variable
, ch.efficiency.sources
### internal variable
, ch.isotope.sigma
### internal variable
, tit
### internal variable
, analysis.name
### internal variable
, filename.prefix
### internal variable
, ch.plot.filename
### internal variable
, plot.format.list
### internal variable
, ch.matrix.sw = 0
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # when ch.matrix.sw = 1, then the matrix plot is performed, otherwise a bunch of individual isotope pairs.

  ### COPY/PASTE FOR DEVEL
  # ch.isotopes.sources      = isotopes.sources      ;
  # ch.isotopes.mixtures     = ch.isotopes.mixtures  ;
  # ch.isotope.mvn.sample    = isotope.mvn.sample    ;
  # ch.concentration.sources = concentration.sources ;
  # ch.efficiency.sources    = efficiency.sources    ;
  # ch.plot.filename         = plot.filename         ;

  #library("chplot"); # will load "ellipse", "KernSmooth", "lattice"
  # sources
  chull.sources                         = matrix(ch.isotopes.sources,      nrow = n.sources);
  chull.concentration.sources           = matrix(ch.concentration.sources, nrow = n.sources);
  chull.efficiency.sources              = matrix(ch.efficiency.sources,    nrow = n.sources);
  colnames(chull.sources)               = names.isotopes; # label columns
  rownames(chull.sources)               = names.sources;  # label rows
  colnames(chull.concentration.sources) = names.isotopes; # label columns
  rownames(chull.concentration.sources) = names.sources;  # label rows
  colnames(chull.efficiency.sources)    = names.isotopes; # label columns
  rownames(chull.efficiency.sources)    = names.sources;  # label rows
  # mixture
  chull.mixture                         = matrix(ch.isotopes.mixtures,      ncol = n.isotopes);
  colnames(chull.mixture)               = names.isotopes; # label columns
  rownames(chull.mixture)               = ch.names.mixtures; # label rows

  color.hull         = "brown";
  color.curved.hull  = "blue";
  color.curved.pairs = "gold1";
  color.sources      = "red";
  color.mixture      = "green";
  color.isotope.mvn  = "orange";
  color.sigma.bars   = "violetred";

  if (ch.n.mixtures == 1) { title.mixture = ch.names.mixtures; };  # if only one mixture, then assign the name of that mixture to the title

  for (i.plot in plot.format.list)
  {
    if (ch.matrix.sw == 1) { # matrix plot
      s.plot.settings.begin.end(filename.prefix, ch.plot.filename, plot.mode = "begin", plot.format = i.plot);
      par(mfrow = c(n.isotopes,n.isotopes));
    } # if ch.matrix.sw=1

    for (i.c in seq(1, n.isotopes)) {
      for (i.r in seq(1, n.isotopes)) {

        # only do upper triangle of matrix when not matrix plot
        if ((ch.matrix.sw == 0) && ((i.r == n.isotopes) || (i.c <= i.r))) { next; }; # if ch.matrix.sw=0

        if (ch.matrix.sw == 0) { # individual plot
          plot.filename = paste(ch.plot.filename, "_", names.isotopes[i.r], "_", names.isotopes[i.c], sep="");
          plot.filename = filename.clean(plot.filename);
          s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "begin", plot.format = i.plot);

          par(mfrow=c(1, 1));
        } # if ch.matrix.sw=0

        # set up axes by finding extremes and expanding them a little for the text labels
        x.min = min(c(chull.sources[, i.r], chull.mixture[, i.r]));
        x.max = max(c(chull.sources[, i.r], chull.mixture[, i.r]));
        y.min = min(c(chull.sources[, i.c], chull.mixture[, i.c]));
        y.max = max(c(chull.sources[, i.c], chull.mixture[, i.c]));
        x.len = x.max-x.min;
        y.len = y.max-y.min;
        x.min = x.min-x.len*0.15; #*0.1;
        x.max = x.max+x.len*0.30; #*0.3;
        y.min = y.min-y.len*0.15; #*0.1;
        y.max = y.max+y.len*0.10; #*0.1;

        plot(NULL, xlab = names.isotopes[i.r], ylab = names.isotopes[i.c], xlim = c(x.min, x.max), ylim = c(y.min, y.max)); # set up plot area and labels
        axis(3); axis(4); # add axis labels to top and right sides

        # plot clouds of isotope.mvn.sample points and std dev bars first
        #for (i.mix.source in seq(1, n.mixtures+n.sources)) {
        #seq.i.mix.source = seq(1, ch.n.mixtures+n.sources);
        for (i.mix.source in seq(1, ch.n.mixtures+n.sources)) {
          # BUG: if there is more than one mixture then not all the sources are plotted and unwanted mixtures are plotted
          # index for groups
          group.index = (i.mix.source-1)*n.isotopes;
          i.r.index = i.r + group.index;
          i.c.index = i.c + group.index;

          # clouds of isotope.mvn.sample points
          #if (i.mix.source == 0) { color.isotope.mvn = color.mixture } else { color.isotope.mvn = color.sources };
          points(x = ch.isotope.mvn.sample[, i.r.index], y = ch.isotope.mvn.sample[, i.c.index], pch = 20, col = color.isotope.mvn, cex = .01);

          # STD DEV BARS
          # coordinates for std dev bars
          arrow.len.denom = 600;  arrow.len.x=x.len/arrow.len.denom;  arrow.len.y=y.len/arrow.len.denom; # lengths of cross
          # horizontal std dev
          if (ch.isotope.sigma[, i.r.index] != 0) {
            y0 = ch.isotope.mvn.sample[, i.c.index];
            x1 = ch.isotope.mvn.sample[, i.r.index] - ch.isotope.sigma[, i.r.index];
            x2 = ch.isotope.mvn.sample[, i.r.index] + ch.isotope.sigma[, i.r.index];
            arrows(x1, y0, x2, y0, length=arrow.len.x, angle=90, code=3, col=color.sigma.bars);
          } # end if ch.isotope.sigma
          # vertical std dev
          if (ch.isotope.sigma[, i.c.index] != 0) {
            x0 = ch.isotope.mvn.sample[, i.r.index];
            y1 = ch.isotope.mvn.sample[, i.c.index] - ch.isotope.sigma[, i.c.index];
            y2 = ch.isotope.mvn.sample[, i.c.index] + ch.isotope.sigma[, i.c.index];
            arrows(x0, y1, x0, y2, length=arrow.len.y, angle=90, code=3, col=color.sigma.bars);
          } # end if ch.isotope.sigma

        } # end for i.mix.source

        if (i.r != i.c) { # non diagonal plots
          # convex hull
          chull.ind.temp = chull (x = chull.sources[, i.r], y = chull.sources[, i.c]);
          chull.ind = c(chull.ind.temp, chull.ind.temp[1]);           # determine indicies of points that make convex hull
          lines(x = chull.sources[chull.ind, i.r], y = chull.sources[chull.ind, i.c], lwd=.5, col = color.hull); # lines making convex hull

          # curved efficiency-concentration between all pairs of points
          curve.spacing = 0.01;
          pairs.ind = t(combn(seq(1,n.sources),2)); # all pairs of indicies in rows
          for (i.line.pairs in seq(1, choose(n.sources,2)) ) {
            line.pair.ind = pairs.ind[i.line.pairs, ]; # indicies of pairs of points

            p.B = matrix(c(seq(0, 1, curve.spacing), seq(1, 0, -curve.spacing)), ncol = 2);  # Biomass proportions
            n.points.curve = dim(p.B)[1];

            # concentrations and efficiencies
            conc.r = chull.concentration.sources[line.pair.ind ,i.r]; # concentration for sources 1 and 2, isotope r
            conc.c = chull.concentration.sources[line.pair.ind ,i.c]; # concentration for sources 1 and 2, isotope c
            effi.r = chull.efficiency.sources[line.pair.ind ,i.r];    # efficiency for sources 1 and 2, isotope r
            effi.c = chull.efficiency.sources[line.pair.ind ,i.c];    # efficiency for sources 1 and 2, isotope c

            # proportions of each isotope for each source
            curved.points = matrix(0, nrow = n.points.curve, ncol = dim(p.B)[2]); # points
            for (i.points.curve in seq(1, n.points.curve)) {
              p.r = p.B[i.points.curve,]*conc.r*effi.r/sum(p.B[i.points.curve,]*conc.r*effi.r); # proportion isotope r
              p.c = p.B[i.points.curve,]*conc.c*effi.c/sum(p.B[i.points.curve,]*conc.c*effi.c); # proportion isotope c
              curved.points[i.points.curve,] = c(sum(chull.sources[line.pair.ind, i.r] * p.r),
                                                 sum(chull.sources[line.pair.ind, i.c] * p.c));
            } # end for i.points.curve

            lines (x = curved.points[, 1], y = curved.points[, 2], lwd=.2, col = color.curved.pairs); # lines making convex hull

          } # end for i.line.pairs curved all pairs

          # curved efficiency-concentration convex hull
          # for each pair of exterior convex hull points in pairs
          curve.spacing = 0.01;
          for (i.hull.pairs in seq(1, length(chull.ind)-1) ) {
            hull.pair.ind = chull.ind[c(i.hull.pairs, i.hull.pairs+1)]; # indicies of pairs of points

            p.B = matrix(c(seq(0, 1, curve.spacing), seq(1, 0, -curve.spacing)), ncol = 2);  # Biomass proportions
            n.points.curve = dim(p.B)[1];

            # concentrations and efficiencies
            conc.r = chull.concentration.sources[hull.pair.ind ,i.r]; # concentration for sources 1 and 2, isotope r
            conc.c = chull.concentration.sources[hull.pair.ind ,i.c]; # concentration for sources 1 and 2, isotope c
            effi.r = chull.efficiency.sources[hull.pair.ind ,i.r];    # efficiency for sources 1 and 2, isotope r
            effi.c = chull.efficiency.sources[hull.pair.ind ,i.c];    # efficiency for sources 1 and 2, isotope c

            # proportions of each isotope for each source
            curved.points = matrix(0, nrow = n.points.curve, ncol = dim(p.B)[2]); # points
            for (i.points.curve in seq(1, n.points.curve)) {
              p.r = p.B[i.points.curve,]*conc.r*effi.r/sum(p.B[i.points.curve,]*conc.r*effi.r); # proportion isotope r
              p.c = p.B[i.points.curve,]*conc.c*effi.c/sum(p.B[i.points.curve,]*conc.c*effi.c); # proportion isotope c
              curved.points[i.points.curve,] = c(sum(chull.sources[hull.pair.ind, i.r] * p.r),
                                                 sum(chull.sources[hull.pair.ind, i.c] * p.c));
            } # end for i.points.curve

            lines (x = curved.points[, 1], y = curved.points[, 2], lwd=2, col = color.curved.hull); # lines making convex hull

          } # end for i.hull.pairs curved convex hull
        } # end if i.r != i.c

        if (ch.matrix.sw == 1) { # matrix plot
          if (i.r == i.c) { # diagonal plots
            # convex hull
            chull.ind.temp = chull (x = chull.sources[, i.r], y = chull.sources[, i.c]);
            chull.ind = c(chull.ind.temp, chull.ind.temp[1]);           # determine indicies of points that make convex hull
            lines (x = chull.sources[chull.ind, i.r], y = chull.sources[chull.ind, i.c], lwd=.5, col = color.hull); # lines making convex hull

            # curved efficiency-concentration between all pairs of points
            curve.spacing = 0.01;
            pairs.ind = t(combn(seq(1,n.sources),2)); # all pairs of indicies in rows
            for (i.line.pairs in seq(1, choose(n.sources,2)) ) {
              line.pair.ind = pairs.ind[i.line.pairs, ]; # indicies of pairs of points

              p.B = matrix(c(seq(0, 1, curve.spacing), seq(1, 0, -curve.spacing)), ncol = 2);  # Biomass proportions
              n.points.curve = dim(p.B)[1];

              # concentrations and efficiencies
              conc.r = chull.concentration.sources[line.pair.ind ,i.r]; # concentration for sources 1 and 2, isotope r
              conc.c = chull.concentration.sources[line.pair.ind ,i.c]; # concentration for sources 1 and 2, isotope c
              effi.r = chull.efficiency.sources[line.pair.ind ,i.r];    # efficiency for sources 1 and 2, isotope r
              effi.c = chull.efficiency.sources[line.pair.ind ,i.c];    # efficiency for sources 1 and 2, isotope c

              # proportions of each isotope for each source
              curved.points = matrix(0, nrow = n.points.curve, ncol = dim(p.B)[2]); # points
              for (i.points.curve in seq(1, n.points.curve)) {
                p.r = p.B[i.points.curve,]*conc.r*effi.r/sum(p.B[i.points.curve,]*conc.r*effi.r); # proportion isotope r
                p.c = p.B[i.points.curve,]*conc.c*effi.c/sum(p.B[i.points.curve,]*conc.c*effi.c); # proportion isotope c
                curved.points[i.points.curve,] = c(sum(chull.sources[line.pair.ind, i.r] * p.r),
                                                   sum(chull.sources[line.pair.ind, i.c] * p.c));
              } # end for i.points.curve

              lines (x = curved.points[, 1], y = curved.points[, 2], lwd=.2, col = color.curved.pairs); # lines making convex hull

            } # end for i.line.pairs curved all pairs

          } # end if i.r == i.c
        } # end if ch.matrix.sw=1

        # points and labels at end to display over lines in body of plot
        points(x = chull.sources[, i.r], y = chull.sources[, i.c], pch = 21, col=color.sources, lwd = 3, cex = 1.2); # plot source points
        points(x = chull.mixture[, i.r], y = chull.mixture[, i.c], pch = 22, col=color.mixture, lwd = 3, cex = 1.2); # mixture point
        text  (x = chull.sources[, i.r], y = chull.sources[, i.c], labels =    names.sources,  adj=c(0,.5), pos = 4, cex = 1  ); # source names of points
        text  (x = chull.mixture[, i.r], y = chull.mixture[, i.c], labels = ch.names.mixtures, adj=c(0,.5), pos = 4, cex = 1.2); # mixture name

        ## points and labels on axes for 1-D view.
        # on horizontal axis
        points(x = chull.sources[, i.r], y = rep(y.min, n.sources),    pch = 21, col=color.sources, lwd = 1, cex = 1); # plot source points
        points(x = chull.mixture[, i.r], y = rep(y.min,ch.n.mixtures), pch = 22, col=color.mixture, lwd = 1, cex = 1); # mixture point
        text  (x = chull.sources[, i.r], y = rep(y.min, n.sources),    labels = substr(names.sources,1,1),     adj=c(0,0), pos = 3, cex = 0.8  ); # source names of points
        text  (x = chull.mixture[, i.r], y = rep(y.min,ch.n.mixtures), labels = substr(ch.names.mixtures,1,1), adj=c(0,0), pos = 3, cex = 1.0); # mixture name
        # on vertical axis
        points(x = rep(x.min, n.sources),    y = chull.sources[, i.c], pch = 21, col=color.sources, lwd = 1, cex = 1); # plot source points
        points(x = rep(x.min,ch.n.mixtures), y = chull.mixture[, i.c], pch = 22, col=color.mixture, lwd = 1, cex = 1); # mixture point
        text  (x = rep(x.min, n.sources),    y = chull.sources[, i.c], labels = substr(   names.sources,1,1),  adj=c(0,0), pos = 4, cex = 0.8  ); # source names of points
        text  (x = rep(x.min,ch.n.mixtures), y = chull.mixture[, i.c], labels = substr(ch.names.mixtures,1,1), adj=c(0,0), pos = 4, cex = 1.0); # mixture name


        if (ch.matrix.sw == 0) { # individual plot
          # outside titles
          s.plot.convex.hull.titles(n.sources, n.isotopes, names.isotopes, analysis.name, title.mixture, tit);
          #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
          s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
        } # end if ch.matrix.sw=0

      } # end for i.r
    } # end for i.c

    if (ch.matrix.sw == 1) { # matrix plot
      # outside titles
      s.plot.convex.hull.titles(n.sources, n.isotopes, names.isotopes, analysis.name, title.mixture, tit);
      #s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "print");
      s.plot.settings.begin.end(filename.prefix, plot.filename, plot.mode = "end");
    } # end if ch.matrix.sw=1

  } # plotting loop

  ### internal variable
} # s.plot.convex.hull.matrix()

# outside titles Convex Hull matrix plot
s.plot.convex.hull.titles = function (n.sources, n.isotopes, names.isotopes, analysis.name, title.mixture, tit)
{
    # box("outer", col="black")  # draw a box around the entire figure
    mtext(paste("Isotopic Mixing Convex Hull"), side=1, line=1, outer=TRUE);                                                # bottom
    mtext(paste("SISUS: Stable Isotope Sourcing using Sampling"), side=1, line=3, outer=TRUE);                              # bottom (line specifies)
    mtext(paste(n.sources, " Sources"), side=2, line=1, outer=TRUE);                                                        # left
    mtext(paste(paste(title.mixture, ": ", tit)), cex = 1.2, side=3, line=1.5, outer=TRUE);                                   # top
    mtext(paste(analysis.name), cex = 0.8, side=3, line=0.25, outer=TRUE);                                                   # top
    mtext(paste(n.isotopes, "Isotopes: ", paste(sprintf("%s",names.isotopes),collapse=", ")), side=4, line=1, outer=TRUE);  # right
    # mtext(expression(paste("Isotopes: ", delta^{13}, "C")), side=4, line=-1, outer=TRUE);  # right
        # having trouble getting the n.isotopes in the front of the line, also trouble feeding the delta^13 as a variable.

  ### internal variable
}
