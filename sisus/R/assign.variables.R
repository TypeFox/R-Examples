assign.variables <-
function# Assign variables from input file to variable names
### internal function for sisus
(DATA
### internal variable
, first=0
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

#
#   When making changes to the workbook, need to:
#     1. get.data.R define sheet names
#     2. assign.variables.R define variables
#     3. sisus.R assign variables
#     4. write.input.R define sheet names
#

    # worksheets from workbook
    SISUS.PARAMETERS                    = DATA$SISUS.PARAMETERS               ;
  if (first != 1) { # don't do if only read first worksheet
    MIXTURES.ISOTOPE.RATIOS             = DATA$MIXTURES.ISOTOPE.RATIOS        ;
    SOURCES.ISOTOPE.RATIOS              = DATA$SOURCES.ISOTOPE.RATIOS         ;
    SOURCES.DISCRIMINATION              = DATA$SOURCES.DISCRIMINATION         ;
    SOURCES.CONCENTRATION               = DATA$SOURCES.CONCENTRATION          ;
    SOURCES.EFFICIENCY                  = DATA$SOURCES.EFFICIENCY             ;
    LINEAR.CONSTRAINTS                  = DATA$LINEAR.CONSTRAINTS             ;
    # prior removed at v0.09 until these ideas are developed 5/22/2008 3:42PM
    PRIOR.DIRICHLET                     = DATA$PRIOR.DIRICHLET                ;
    MIXTURES.ISOTOPE.STDDEV             = DATA$MIXTURES.ISOTOPE.STDDEV        ;
    SOURCES.ISOTOPE.STDDEV              = DATA$SOURCES.ISOTOPE.STDDEV         ;
    SOURCES.BIOMASS.PER.INDIVIDUAL      = DATA$SOURCES.BIOMASS.PER.INDIVIDUAL ;
    SOURCES.NUMBER.OF.INDIVIDUALS       = DATA$SOURCES.NUMBER.OF.INDIVIDUALS  ;
    #MIXTURES.CONCENTRATION.STDDEV       = DATA$MIXTURES.CONCENTRATION.STDDEV  ;
    #SOURCES.CONCENTRATION.STDDEV        = DATA$SOURCES.CONCENTRATION.STDDEV   ;
    #MIXTURES.EFFICIENCY.STDDEV          = DATA$MIXTURES.EFFICIENCY.STDDEV     ;
    #SOURCES.EFFICIENCY.STDDEV           = DATA$SOURCES.EFFICIENCY.STDDEV      ;
  } # if first

  ## SISUS.PARAMETERS       ######################################################
  # parameter settings
  param.values = as.vector(type.convert( SISUS.PARAMETERS[seq(1,dim(SISUS.PARAMETERS)[1]), 2] ));

    ### Column 1 names  ============================  Row  =====================
    ## SISUS Parameters
    # Analysis Labels
    row.workbook =  2;
    row.workbook = row.workbook + 1;
    analysis.name                   =               param.values[row.workbook];  row.workbook = row.workbook + 1;
    filename.prefix                 =               param.values[row.workbook];  row.workbook = row.workbook + 1;
    mixtures.group.name             =               param.values[row.workbook];  row.workbook = row.workbook + 1;
      if (identical(filename.prefix, "")) { filename.prefix =  "SISUS"; } # if no filename.prefix supplied

    # Sources and Isotopes
    row.workbook = row.workbook + 1;
    n.mixtures                      = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    n.sources                       = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    n.isotopes                      = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;

    # Sampling
    row.workbook = row.workbook + 1;
    M                               = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    ##n.samples.isotope.mvn           = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
      n.samples.isotope.mvn = 0;
    skip                            = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    burnin                          = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    seed                            = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;

    SW = new.env();  # create a list for Mixing Model and Output Settings switches

    # Mixing Model
    row.workbook = row.workbook + 1;
    simplex.include.sw                = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    discrimination.include.sw         = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    concentration.include.sw          = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    assimeffic.include.sw             = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    lc.include.sw                     = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    # prior removed at v0.09 until these ideas are developed 5/22/2008 3:42PM
    #prior.include.sw                  = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    isotope.stddev.include.sw         = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    biomass.per.individual.include.sw = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    ##number.of.individuals.include.sw  = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    #concentration.stddev.include.sw   = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    #efficiency.stddev.include.sw      = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;

        SW$simplex.include.sw                     = simplex.include.sw;
        SW$discrimination.include.sw              = discrimination.include.sw;
        SW$concentration.include.sw               = concentration.include.sw;
        SW$assimeffic.include.sw                  = assimeffic.include.sw;
        SW$lc.include.sw                          = lc.include.sw;
        # prior removed at v0.09 until these ideas are developed 5/22/2008 3:42PM
          prior.include.sw                        = 0;
        SW$prior.include.sw                       = prior.include.sw;
        SW$isotope.stddev.include.sw              = isotope.stddev.include.sw;
        SW$biomass.per.individual.include.sw      = biomass.per.individual.include.sw;
        ##SW$number.of.individuals.include.sw       = number.of.individuals.include.sw;
        #SW$concentration.stddev.include.sw        = concentration.stddev.include.sw;
        #SW$efficiency.stddev.include.sw           = efficiency.stddev.include.sw;

    # Output Settings
    row.workbook = row.workbook + 1;
    plot.convex.hull.include.sw          = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.marginal.histogram.include.sw   = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.scatterplot.matrix.include.sw   = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    numerical.summaries.include.sw       = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    solution.polytope.vertices.sw        = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    samples.include.sw                   = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    mcmc.diagnostics.sw                  = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    time.series.sw                       = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    hist.bins                            = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    hist.density.sw                      = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;

        SW$plot.convex.hull.include.sw          = plot.convex.hull.include.sw       ;
        SW$plot.marginal.histogram.include.sw   = plot.marginal.histogram.include.sw;
        SW$plot.scatterplot.matrix.include.sw   = plot.scatterplot.matrix.include.sw;
        SW$numerical.summaries.include.sw       = numerical.summaries.include.sw    ;
        SW$solution.polytope.vertices.sw        = solution.polytope.vertices.sw     ;
        SW$samples.include.sw                   = samples.include.sw                ;
        SW$mcmc.diagnostics.sw                  = mcmc.diagnostics.sw               ;
        SW$time.series.sw                       = time.series.sw                    ;
        SW$hist.density.sw                      = hist.density.sw                   ;

    # Plot Formats
    row.workbook = row.workbook + 1;
    plot.format.png                 = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.format.eps                 = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.format.pdf                 = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.format.bmp                 = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.format.jpeg                = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    #plot.format.gif                 = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;
    plot.format.tiff                = type.convert( param.values[row.workbook]); row.workbook = row.workbook + 1;

      # vector of plot types desired
      plot.format.list = NULL;
      #if (plot.format.png  != 0) { plot.format.list = c(plot.format.list, 1); }
                                   plot.format.list = c(plot.format.list, 1); # always produce png since bmp, jpeg, tiff requires under unix
      if (plot.format.eps  != 0) { plot.format.list = c(plot.format.list, 2); }
      if (plot.format.pdf  != 0) { plot.format.list = c(plot.format.list, 3); }
      if (plot.format.bmp  != 0) { plot.format.list = c(plot.format.list, 4); }
      if (plot.format.jpeg != 0) { plot.format.list = c(plot.format.list, 5); }
    #  if (plot.format.gif  != 0) { plot.format.list = c(plot.format.list, 6); }
      if (plot.format.tiff != 0) { plot.format.list = c(plot.format.list, 6); }

  ########################################
  if (first == 1) { # return switches if first worksheet
    return(SW);
  }

  ## MIXTURES.ISOTOPE.RATIOS        ######################################################
  # names of source, isotopes, and sources
  names.isotopes = MIXTURES.ISOTOPE.RATIOS[1,seq(2,1+n.isotopes)];
  names.mixtures = MIXTURES.ISOTOPE.RATIOS[seq(2,1+n.mixtures),1];

  # observed number of mixtures and isotopes in workbook
  #n.mixtures.worksheet = dim(MIXTURES.ISOTOPE.RATIOS)[1]-1;
  n.mixtures.worksheet = sum((MIXTURES.ISOTOPE.RATIOS[,1] != "NA"), na.rm=TRUE)-1; # count the number of non NA fields and subtract header
  if (n.mixtures != n.mixtures.worksheet) {
          p.o = paste("\n"); write.out(p.o);
    p.o = paste("Information:  Using", n.mixtures, "mixtures, but there appears to be", n.mixtures.worksheet, "in worksheet.", "\n"); write.out(p.o);
  }
  #n.isotopes.worksheet = dim(MIXTURES.ISOTOPE.RATIOS)[2]-1;
  n.isotopes.worksheet = sum((MIXTURES.ISOTOPE.RATIOS[1,] != "NA"), na.rm=TRUE)-1; # count the number of non NA fields and subtract header
  if (n.isotopes != n.isotopes.worksheet) {
          p.o = paste("\n"); write.out(p.o);
    p.o = paste("Information:  Using", n.isotopes, "isotopes, but there appears to be", n.isotopes.worksheet, "in worksheet.", "\n"); write.out(p.o);
  }

    # an extra space is added to end of last isotope name (unknown why -- artifact of reading Excel workbook)
    # remove that last space if it exists
    if (substr(names.isotopes[n.isotopes], nchar(names.isotopes[n.isotopes]), nchar(names.isotopes[n.isotopes])) == " ") {
      names.isotopes[n.isotopes] = substr(names.isotopes[n.isotopes],1,nchar(names.isotopes[n.isotopes])-1);
    }

  # initial isotope values
  isotopes.init.mixtures = type.convert( MIXTURES.ISOTOPE.RATIOS[seq(2,1+n.mixtures),seq(2,1+n.isotopes)] );

  ## SOURCES.ISOTOPE.RATIOS        ######################################################
  # names of source, isotopes, and sources
  names.sources  = SOURCES.ISOTOPE.RATIOS[seq(2,1+n.sources),1];

  #n.sources.worksheet = dim(SOURCES.ISOTOPE.RATIOS)[2]-1;
  n.sources.worksheet = sum((SOURCES.ISOTOPE.RATIOS[,1] != "NA"), na.rm=TRUE)-1; # count the number of non NA fields and subtract header
  if (n.sources != n.sources.worksheet) {
          p.o = paste("\n"); write.out(p.o);
    p.o = paste("Information:  Using", n.sources, "sources,  but there appears to be", n.sources.worksheet, "in worksheet.", "\n"); write.out(p.o);
  }

  # initial isotope values
  isotopes.init.sources = type.convert( SOURCES.ISOTOPE.RATIOS[seq(2,1+n.sources),seq(2,1+n.isotopes)] );

  ## SOURCES.DISCRIMINATION       ######################################################
  # discrimination values
  if (SW$discrimination.include.sw == 1) {
    discrimination.sources = type.convert( SOURCES.DISCRIMINATION[seq(2,1+n.sources),seq(2,1+n.isotopes)] );
  } else {
    discrimination.sources = matrix(data = 0, nrow=n.sources, ncol=n.isotopes);  # no concentration
  }

    isotopes.mixtures = isotopes.init.mixtures;
  # discrimination corrected isotope signatures
    isotopes.sources = isotopes.init.sources + discrimination.sources;

  ## SOURCES.CONCENTRATION        ######################################################
  # concentration values
  if (SW$concentration.include.sw == 1) {
    concentration.sources = type.convert( SOURCES.CONCENTRATION[seq(2,1+n.sources),seq(2,1+n.isotopes)] );
  } else {
    concentration.sources = matrix(data = 1, nrow=n.sources, ncol=n.isotopes);  # no concentration
  }

  ## SOURCES.EFFICIENCY           ######################################################
  # assimilation efficiency values
  if (SW$assimeffic.include.sw == 1) {
    efficiency.sources = type.convert( SOURCES.EFFICIENCY[seq(2,1+n.sources),seq(2,1+n.isotopes)] );
  } else {
    efficiency.sources = matrix(data = 1, nrow=n.sources, ncol=n.isotopes);  # no efficiency
  }

  ## LINEAR.CONSTRAINTS          ######################################################
  # additional linear constraints
  if (SW$lc.include.sw == 1) {
    constraints.cols    = dim(LINEAR.CONSTRAINTS)[2];
    constraints.type    = trim(as.vector( LINEAR.CONSTRAINTS[1,seq(2,constraints.cols)] ));
    constraints.RHS.b   = type.convert( LINEAR.CONSTRAINTS[2,seq(2,constraints.cols)] );
    constraints.sources = type.convert( LINEAR.CONSTRAINTS[seq(3,2+n.sources),seq(2,constraints.cols)] );
  }

    # define additional linear constraints
    lc = additional.linear.constraints(constraints.cols, constraints.type, constraints.RHS.b, constraints.sources, n.sources, SW$lc.include.sw);

  ## PRIOR.DIRICHLET               ######################################################
  # priors
  if (SW$prior.include.sw == 1) { # if no prior, then specify as uniform
    #priors.cols      = dim(PRIOR.DIRICHLET)[2];
    priors.cols      = 2;  # hardset to 2 since data only in first two columns.  It was causing an error by using 12 columns of NA's.
    priors.type      = trim(as.vector( PRIOR.DIRICHLET[1,seq(2,priors.cols)] ));
    priors.precision = type.convert(   PRIOR.DIRICHLET[2,seq(2,priors.cols)] );
    priors.sources   = type.convert(   PRIOR.DIRICHLET[seq(3,2+n.sources),seq(2,priors.cols)] );
  } else {
    priors.cols      = 2;  # hardset to 2 since data only in first two columns.  It was causing an error by using 12 columns of NA's.
    priors.type      = "p proportions";
    priors.precision = 1;
    priors.sources = rep(1/n.sources, times=n.sources);
  }

  # isotope standard deviations
  ## MIXTURES.ISOTOPE.STDDEV       ######################################################
  if (SW$isotope.stddev.include.sw == 1) {
    isotope.sigma.mixtures = type.convert( MIXTURES.ISOTOPE.STDDEV[seq(2,1+n.mixtures),seq(2,1+n.isotopes)] );
  } else {
    isotope.sigma.mixtures = matrix(data = 0, nrow=n.mixtures, ncol=n.isotopes);  # no variation
  }

  ## SOURCES.ISOTOPE.STDDEV       ######################################################
  if (SW$isotope.stddev.include.sw == 1) {
    isotope.sigma.sources = type.convert( SOURCES.ISOTOPE.STDDEV[seq(2,1+n.sources),seq(2,1+n.isotopes)] );
  } else {
    isotope.sigma.sources = matrix(data = 0, nrow=n.sources, ncol=n.isotopes);  # no variation
  }

    # construct mean vector and covariance matrix for isotopes
    isotope.mean   = t(as.matrix(c(t(isotopes.mixtures), t(isotopes.sources))));
    isotope.sigma = diag(c(t(isotope.sigma.mixtures), t(isotope.sigma.sources))); # diagonal covariance matrix


  # Proportions of Sources Consumed model
  ## SOURCES.BIOMASS.PER.INDIVIDUAL ###################################################
  if (SW$biomass.per.individual.include.sw == 1) {
    biomass.per.individual.sources = matrix(data = rep( type.convert( SOURCES.BIOMASS.PER.INDIVIDUAL[seq(2,1+n.sources),2] ) , n.isotopes), ncol=n.isotopes, byrow = FALSE);
    number.of.individuals.sources  = matrix(data = rep( type.convert( SOURCES.BIOMASS.PER.INDIVIDUAL[seq(2,1+n.sources),3] ) , n.isotopes), ncol=n.isotopes, byrow = FALSE);
  } else {
    biomass.per.individual.sources = matrix(data = 1, nrow=n.sources, ncol=n.isotopes);  # not used
    number.of.individuals.sources  = matrix(data = 1, nrow=n.sources, ncol=n.isotopes);  # not used
  }


 ## # concentration standard deviations
 ## ## MIXTURES.CONCENTRATION.STDDEV       ######################################################
 ## if (SW$concentration.stddev.include.sw == 1) {
 ##   concentration.sigma.mixtures = type.convert( MIXTURES.CONCENTRATION.STDDEV[seq(2,1+n.mixtures),seq(2,1+n.isotopes)] );
 ## } else {
 ##   concentration.sigma.mixtures = matrix(data = 0, nrow=n.mixtures, ncol=n.isotopes);  # no variation
 ## }
 ##
 ## ## SOURCES.CONCENTRATION.STDDEV       ######################################################
 ## if (SW$concentration.stddev.include.sw == 1) {
 ##   concentration.sigma.sources = type.convert( SOURCES.CONCENTRATION.STDDEV[seq(2,1+n.sources),seq(2,1+n.isotopes)] );
 ## } else {
 ##   concentration.sigma.sources = matrix(data = 0, nrow=n.sources, ncol=n.isotopes);  # no variation
 ## }
 ##
 ##   # construct mean vector and covariance matrix for isotopes
 ##   concentration.mean   = t(as.matrix(c(t(isotopes.mixtures), t(isotopes.sources))));
 ##   concentration.sigma = diag(c(t(isotope.sigma.mixtures), t(isotope.sigma.sources))); # diagonal covariance matrix
 ##
 ##
 ## # efficiency standard deviations
 ## ## MIXTURES.EFFICIENCY.STDDEV       ######################################################
 ## if (SW$efficiency.stddev.include.sw == 1) {
 ##   efficiency.sigma.mixtures = type.convert( MIXTURES.EFFICIENCY.STDDEV[seq(2,1+n.mixtures),seq(2,1+n.isotopes)] );
 ## } else {
 ##   efficiency.sigma.mixtures = matrix(data = 0, nrow=n.mixtures, ncol=n.isotopes);  # no variation
 ## }
 ##
 ## ## SOURCES.EFFICIENCY.STDDEV       ######################################################
 ## if (SW$efficiency.stddev.include.sw == 1) {
 ##   efficiency.sigma.sources = type.convert( SOURCES.EFFICIENCY.STDDEV[seq(2,1+n.sources),seq(2,1+n.isotopes)] );
 ## } else {
 ##   efficiency.sigma.sources = matrix(data = 0, nrow=n.sources, ncol=n.isotopes);  # no variation
 ## }
 ##
 ##   # construct mean vector and covariance matrix for isotopes
 ##   efficiency.mean   = t(as.matrix(c(t(isotopes.mixtures), t(isotopes.sources))));
 ##   efficiency.sigma = diag(c(t(isotope.sigma.mixtures), t(isotope.sigma.sources))); # diagonal covariance matrix


  ##############################################################################
  VARIABLES = new.env();  # create a list to return with all parameters to run SISUS

  # Analysis Labels
  VARIABLES$analysis.name                          = analysis.name;
  VARIABLES$filename.prefix                        = filename.prefix;
  VARIABLES$mixtures.group.name                    = mixtures.group.name;

  # Data Dimensions
  VARIABLES$n.mixtures                             = n.mixtures;
  VARIABLES$n.sources                              = n.sources;
  VARIABLES$n.isotopes                             = n.isotopes;

  # Sampling
  VARIABLES$M                                      = M;
  VARIABLES$n.samples.isotope.mvn                  = n.samples.isotope.mvn;
  VARIABLES$skip                                   = skip;
  VARIABLES$burnin                                 = burnin;
  VARIABLES$seed                                   = seed;

  # Mixing Model switches
  VARIABLES$SW                                     = SW;

  # Plot Settings
  VARIABLES$hist.bins                              = hist.bins;

  # Plot Formats
  VARIABLES$plot.format.list                       = plot.format.list;

  # names of source, isotopes, and sources
  VARIABLES$names.isotopes                         = names.isotopes;
  VARIABLES$names.mixtures                         = names.mixtures;
  VARIABLES$names.sources                          = names.sources;

  # discrimination corrected isotope signatures
  VARIABLES$isotopes.mixtures                      = isotopes.mixtures;
  VARIABLES$isotopes.sources                       = isotopes.sources;

  # concentration values
  VARIABLES$concentration.sources                  = concentration.sources;

  # assimilation efficiency values
  VARIABLES$efficiency.sources                     = efficiency.sources;

  # define additional linear constraints
  VARIABLES$lc                                     = lc;

  # priors
  VARIABLES$priors.cols                            = priors.cols;
  VARIABLES$priors.type                            = priors.type;
  VARIABLES$priors.precision                       = priors.precision;
  VARIABLES$priors.sources                         = priors.sources;

  # construct mean vector and covariance matrix for isotopes
  VARIABLES$isotope.mean                           = isotope.mean;
  VARIABLES$isotope.sigma                          = isotope.sigma;

  # biomass per individual values
  VARIABLES$biomass.per.individual.sources         = biomass.per.individual.sources;

  # number of individuals values
  VARIABLES$number.of.individuals.sources          = number.of.individuals.sources;

 # # construct mean vector and covariance matrix for concentration
 # VARIABLES$concentration.mean                     = concentration.mean;
 # VARIABLES$concentration.sigma                    = concentration.sigma;
 #
 # # construct mean vector and covariance matrix for efficiency
 # VARIABLES$efficiency.mean                        = efficiency.mean;
 # VARIABLES$efficiency.sigma                       = efficiency.sigma;

  return( as.list(VARIABLES) );
  ### internal variable
}
