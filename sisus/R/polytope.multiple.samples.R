polytope.multiple.samples <-
function# Defines linear constraints for each sample from the isotope distribution, passed to sample.from.polytope()
### internal function for sisus
(M
### internal variable
, n.samples.isotope.mvn
### internal variable
, n.isotopes
### internal variable
, n.sources
### internal variable
, isotope.mvn.sample.indy.mixture
### internal variable
, concentration.sources
### internal variable
, efficiency.sources
### internal variable
, biomass.per.individual.sources
### internal variable
, number.of.individuals.sources
### internal variable
, lc
### internal variable
, skip
### internal variable
, burnin
### internal variable
, names.sources
### internal variable
, SW
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  ########################################
  # Draw from polytope
  # for each sample from the distribution of isotopes
  # use those values to create a polytope and draw some samples from that polytope
  # append all the samples
  sam = NULL;
  sol.feasible.count = 0;
  warning.sw = 0; # only give one warning for unfeasible solutions
  n.samples.isotope.mvn.polytope = ceiling(M/n.samples.isotope.mvn); # number to sample from each polytope
        p.o = paste("           Drawing", n.samples.isotope.mvn, "sample(s) of size", n.samples.isotope.mvn.polytope, "for a total of ", M, "samples: "); write.out(p.o);
  for (i.samples.isotope.mvn in seq(1, n.samples.isotope.mvn)) { # number of polytopes to sample from
        if ((i.samples.isotope.mvn %% 20) == 1) { p.o = paste("\n           "); write.out(p.o); } # "%%" is mod, 20 samples to a line
        p.o = paste(i.samples.isotope.mvn, " "); write.out(p.o);
    # for each isotope.mvn.sample.indy.mixture, break it into mixture and sources matricies
    isotopes.all.temp = matrix(isotope.mvn.sample.indy.mixture[i.samples.isotope.mvn,], ncol = n.isotopes, byrow = TRUE);
    isotopes.mixtures.temp = isotopes.all.temp[1,];
    isotopes.sources.temp = isotopes.all.temp[seq(2, n.sources+1),];

    # Model Mass Balance Equation
    # efficiency-concentration-dependent biomass equations,
    #   subtract mixture and multiply by efficiency and concentration: ( dE_{i',j} - dE_{M,j} ) * [E_{i,j}]
    biomass.values.mvn = model.mass.balance.equation(isotopes.mixtures.temp, isotopes.sources.temp, concentration.sources, efficiency.sources, biomass.per.individual.sources, number.of.individuals.sources);

    # set up linear constraints based on isotope sample and additional linear constraints
    Ab = polytope.constraints(lc, n.sources, n.isotopes, biomass.values.mvn, SW$simplex.include.sw); # define polytope linear constraints

    # Write out Ab matricies as equations in log file
    write.Ab(Ab);

    # CORE FUNCTION
    # sample from the polytope
    SAMPLE = sample.from.polytope(Ab, M, skip, burnin, warning.sw, i.samples.isotope.mvn);  # CORE FUNCTION

    if (i.samples.isotope.mvn == 1) {
      ##V.sam      = SAMPLE$V.sam;      # vertices of original polytope
      n.vertices = SAMPLE$n.vertices; # number of vertices of original polytope
    }

    warning.sw = SAMPLE$warning.sw; # update warning switch
    if (SAMPLE$sol.feasible == 1) {
      sam = rbind(sam, SAMPLE$sam);  # append samples
      sol.feasible.count = sol.feasible.count + 1;  # count feasible solutions
    };
  } # end for i.samples.isotope.mvn
        # new line after printing sample numbers
        p.o = paste("\n");  write.out(p.o);

  if (is.null(sam)) {
    M.actual = 0;
  } else {
    M.actual = dim(sam)[1]-n.vertices;    # reset sample size M with actual sample size (ie., set to 1 if unique)
    colnames(sam) = names.sources;        # label columns with source names
  }

  POLYTOPE.SAMPLE = new.env();  # create a list to return
  POLYTOPE.SAMPLE$M.actual           = M.actual          ;
  POLYTOPE.SAMPLE$sam                = sam               ;
  POLYTOPE.SAMPLE$sol.feasible.count = sol.feasible.count;
  POLYTOPE.SAMPLE$n.vertices         = n.vertices        ;
  ##POLYTOPE.SAMPLE$V.sam              = V.sam             ;
  return( as.list(POLYTOPE.SAMPLE) );
  ### internal variable
}
