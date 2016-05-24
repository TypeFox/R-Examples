model.mass.balance.equation <-
function# Mass balance equations for converting isotope values to biomass values
### internal function for sisus
(isotopes.mixtures
### internal variable
, isotopes.sources
### internal variable
, concentration.sources
### internal variable
, efficiency.sources
### internal variable
, biomass.per.individual.sources
### internal variable
, number.of.individuals.sources
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  biomass.values = t(t(isotopes.sources) - isotopes.mixtures) * (biomass.per.individual.sources * number.of.individuals.sources) * concentration.sources * efficiency.sources;
  return(biomass.values);
  ### internal variable
}
