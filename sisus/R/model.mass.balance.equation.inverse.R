model.mass.balance.equation.inverse <-
function# Inverse Mass balance equations for converting biomass sampled values back to isotope values
### internal function for sisus
(M
### internal variable
, n.vertices
### internal variable
, n.isotopes
### internal variable
, n.sources
### internal variable
, p.biomass.sam
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

  ## Back-solve for p.isotopes

  # p.isotopes = matrix(0,n.sources,n.isotopes); # init to 0s
  # for (i in seq(1,n.sources)) {
  #   for (j in seq(1,n.isotopes)) {
  #     p.isotopes[i,j] = p.biomass.sam[i]*concentration.sources[i,j] / sum( p.biomass.sam*concentration.sources[,j] );
  #   }
  # }

  # make them matricies if they are only vectors -- when n.sources=2
  efficiency.sources             = as.matrix(efficiency.sources);
  concentration.sources          = as.matrix(concentration.sources);
  biomass.per.individual.sources = as.matrix(biomass.per.individual.sources);
  number.of.individuals.sources  = as.matrix(number.of.individuals.sources);

  # p.isotopes has M rows, the columns are groups of sources for each element isotope {iso1(s1, .., sn), .., ison(s1, .., sn)}
  p.isotopes = matrix(0, n.vertices+M, n.sources*n.isotopes); # init to 0s
  for (it.M in seq(1, n.vertices+M)) {
    for (j in seq(1, n.isotopes)) {
      for (i in seq(1, n.sources)) {
        ind.col.ij = i+(j-1)*n.sources;
        p.biomass.num =      as.numeric(p.biomass.sam[it.M, i]) * (biomass.per.individual.sources[i, j] * number.of.individuals.sources[i, j]) * concentration.sources[i, j] * efficiency.sources[i, j]  ;
        p.biomass.den = sum( as.numeric(p.biomass.sam[it.M,  ]) * (biomass.per.individual.sources[ , j] * number.of.individuals.sources[ , j]) * concentration.sources[ , j] * efficiency.sources[ , j] );
        p.isotopes[it.M, ind.col.ij] = p.biomass.num / p.biomass.den ;

        #print(c(it.M, j, i, p.biomass.num, p.biomass.den))
      }
    }
  }

  # if a polytope vertex is all 0s (as possible in PSC model), then value will be NaN -- replace with 0s
  p.isotopes[is.nan(p.isotopes)]=0;
  p.isotopes[is.na(p.isotopes)]=0;

  return(p.isotopes);
  ### internal variable
}
