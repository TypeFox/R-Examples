write.model.settings <-
function# Write model settings to the process log
### internal function for sisus
(n.mixtures
### internal variable
, n.sources
### internal variable
, n.isotopes
### internal variable
, analysis.name
### internal variable
, filename.prefix
### internal variable
, mixtures.group.name
### internal variable
, SW
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()


        p.o = paste("\n"); write.out(p.o);
        p.o = paste("### SISUS Model Settings ###", "\n"); write.out(p.o);

                                                    p.o = paste("  analysis.name = ",       analysis.name, "\n");       write.out(p.o);
                                                    p.o = paste("  filename.prefix = ",     filename.prefix, "\n");     write.out(p.o);
                                                    p.o = paste("  mixtures.group.name = ", mixtures.group.name, "\n"); write.out(p.o);
        p.o = paste("\n"); write.out(p.o);
                                                    p.o = paste("  m = ", n.mixtures, " Mixtures", "\n"); write.out(p.o);
                                                    p.o = paste("  n = ", n.sources,  " Sources",  "\n"); write.out(p.o);
                                                    p.o = paste("  k = ", n.isotopes, " Isotopes", "\n"); write.out(p.o);
        p.o = paste("\n"); write.out(p.o);
  if (SW$simplex.include.sw                == 1)  { p.o = paste("  Include SIMPLEX (sum of proportions is 1) constraint -- for Mixing Model", "\n"); }
                                             else { p.o = paste("  NOT including SIMPLEX -- for Proportion of Sources Consumed Mixing Model", "\n"); };  write.out(p.o);
        p.o = paste("\n"); write.out(p.o);
  if (SW$discrimination.include.sw         == 1)  { p.o = paste("  Correcting source isotopic signatures for discrimination", "\n"); }
                                             else { p.o = paste("  NOT correcting source isotopic signatures for discrimination", "\n"); };                write.out(p.o);
  if (SW$concentration.include.sw          == 1)  { p.o = paste("  Using SOURCE CONCENTRATIONS", "\n"); }
                                             else { p.o = paste("  NOT using SOURCE CONCENTRATIONS", "\n"); };                                             write.out(p.o);
  if (SW$assimeffic.include.sw             == 1)  { p.o = paste("  Using ASSIMILATION EFFICIENCIES", "\n"); }
                                             else { p.o = paste("  NOT using ASSIMILATION EFFICIENCIES", "\n"); };                                         write.out(p.o);
  if (SW$lc.include.sw                     == 1)  { p.o = paste("  Including additional LINEAR CONSTRAINTS", "\n"); }
                                             else { p.o = paste("  NOT including additional LINEAR CONSTRAINTS", "\n"); };                                 write.out(p.o);
  if (SW$prior.include.sw                  == 1)  { p.o = paste("  Importance resampling using PRIOR weighting of uniform samples, if specified, in addition to uniform samples from solution polytope", "\n"); }
                                             else { p.o = paste("  NOT using PRIOR, only uniform samples from solution polytope", "\n"); };                write.out(p.o);
  if (SW$isotope.stddev.include.sw         == 1)  { p.o = paste("  Including ISOTOPE STDDEV on Convex Hull plots", "\n"); }
                                             else { p.o = paste("  NOT using ISOTOPE STDDEV", "\n"); };                                                    write.out(p.o);
  if (SW$biomass.per.individual.include.sw == 1 && SW$simplex.include.sw == 0)  { p.o = paste("  Including BIOMASS PER INDIVIDUAL", "\n"); }
                                             else { p.o = paste("  NOT using BIOMASS PER INDIVIDUAL", "\n"); };                                              write.out(p.o);
  if (SW$biomass.per.individual.include.sw == 1 && SW$simplex.include.sw == 0)  { p.o = paste("  Including NUMBER OF INDIVIDUALS", "\n"); }
                                             else { p.o = paste("  NOT using NUMBER OF INDIVIDUALS", "\n"); };                                                 write.out(p.o);
 # if (SW$concentration.stddev.include.sw   == 1)  { p.o = paste("  ** IGNORED, not yet implemented **", "\n"); }
 #                                            else { p.o = paste("  NOT using CONCENTRATION STDDEV", "\n"); };                                              write.out(p.o);
 # if (SW$efficiency.stddev.include.sw      == 1)  { p.o = paste("  ** IGNORED, not yet implemented **", "\n"); }
 #                                            else { p.o = paste("  NOT using EFFICIENCY STDDEV", "\n"); };                                                 write.out(p.o);






        p.o = paste("\n"); write.out(p.o);
        p.o = paste("### SISUS Model Assumptions (under development) ###", "\n"); write.out(p.o);

  if (SW$simplex.include.sw                == 1)  { p.o = paste("  Mixing Model using Simplex", "\n"); }
                                             else { p.o = paste("  Proportions of Sources Consumed Mixing Model", "\n"); };  write.out(p.o);
  if (SW$discrimination.include.sw         == 1)  { p.o = paste("  Discrimination values are accurately specified.", "\n"); }
                                             else { p.o = paste("  No discrimination or isotope values already corrected.", "\n"); };  write.out(p.o);
  if (SW$concentration.include.sw          == 1)  { p.o = paste("  Elemental concentrations are accurately specified.", "\n"); }
                                             else { p.o = paste("  Elemental concentrations are the same for all sources and elements.", "\n"); };                               write.out(p.o);
  if (SW$assimeffic.include.sw             == 1)  { p.o = paste("  Assimilation efficiencies are accurately specified.", "\n"); }
                                             else { p.o = paste("  Assimilation efficiencies are the same for all sources and elements.", "\n"); };                           write.out(p.o);
  if (SW$lc.include.sw                     == 1)  { p.o = paste("  Additional linear constraints are accurately specified.", "\n"); }
                                             else { p.o = NULL; };                   write.out(p.o);
  if (SW$prior.include.sw                  == 1)  { p.o = paste("  Prior (expert) knowledge properly specified as Dirichlet distribution over the simplex.", "\n"); }
                                             else { p.o = paste("  Solutions within solution polytope are equally likely.", "\n"); };  write.out(p.o);
  #if (SW$isotope.stddev.include.sw         == 1)  { p.o = paste("  ** Experimental ** Variance in isotope ratios, following independent normal distributions with specified center and standard deviation.", "\n"); }
  #                                           else { p.o = paste("  No variance in isotope ratios in sources.", "\n"); };                                      write.out(p.o);
  if (SW$biomass.per.individual.include.sw == 1 && SW$simplex.include.sw == 0)  { p.o = paste("  Biomass per individual accurately specified, using Proportions of Sources Consumed Mixing Model", "\n"); }
                                             else { p.o = paste("  NOT using BIOMASS PER INDIVIDUAL", "\n"); };                                              write.out(p.o);
  if (SW$biomass.per.individual.include.sw == 1 && SW$simplex.include.sw == 0)  { p.o = paste("  Number of Individuals in consumable population accurately specified, using Proportions of Sources Consumed Mixing Model", "\n"); }
                                             else { p.o = paste("  NOT using NUMBER OF INDIVIDUALS", "\n"); };                                                 write.out(p.o);
 # if (SW$concentration.stddev.include.sw   == 1)  { p.o = paste("  No variance in elemental concetrations in sources. ** not yet implemented **", "\n"); }
 #                                            else { p.o = paste("  No variance in elemental concetrations in sources.", "\n"); };                                write.out(p.o);
 # if (SW$efficiency.stddev.include.sw      == 1)  { p.o = paste("  No variance in assimilation efficiencies. ** not yet implemented **", "\n"); }
 #                                            else { p.o = paste("  No variance in assimilation efficiencies.", "\n"); };                                   write.out(p.o);

        p.o = paste("\n"); write.out(p.o);
        p.o = paste("### SISUS Output Settings ###", "\n"); write.out(p.o);

        p.o = paste("\n"); write.out(p.o);
  if (SW$time.series.sw                    == 1)  { p.o = paste("  Treating Mixtures as a time series", "\n"); }
                                             else { p.o = paste("  Mixtures are not a time series", "\n"); };                                              write.out(p.o);
  if (SW$time.series.sw                    == 1)  { p.o = paste("  Mixtures as a time series assumes sources are not changing", "\n"); }
                                             else { p.o = NULL; };                                                 write.out(p.o);

        p.o = paste("\n"); write.out(p.o);


  ### internal variable
}
