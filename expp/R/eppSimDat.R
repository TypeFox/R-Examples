
eppSimDat <-function(N = 10, meanClutch = 10, eppRate = 0.10, eppMax = 12, eppMales = 0.35, nLags = 3) {
  
  # breeding
  d = data.frame(x = rnorm(N), y = rnorm(N), id = 1:N, male = paste0("m", 1:N), female = paste0("f", 1:N), stringsAsFactors = F )
  d$clutch = rpois(N, meanClutch)
  d$trait = rnorm(N)
  
  # epp
  d$epy = round(d$clutch * rnorm(N, mean = eppRate, sd = eppRate)    )
  d$epy = ifelse(d$epy > d$clutch, d$clutch, ifelse(d$epy < 0, 0, d$epy) )
  
  d$epmale = FALSE
  d[ sample(1:N, round(N*eppMales) ), "epmale"] = TRUE
  
  # epp pairs
  females = rep(d$female,  times = d$epy)
  males = d[d$epmale, "male"]
  males = sapply(females, function(f)    sample( setdiff(males, d[d$female == f, "male"] ), 1) )
  
  eppPairs = data.frame(female = females, male = males)
  
  e= eppMatrix(eppPairs, pairs = ~ male + female)
  
  # the epp object
  d = SpatialPointsBreeding(d[, c("x", "y", "male", "female", "id", "trait")], id = 'id')
  
  epp(d, DirichletPolygons(d) , e, maxlag = nLags)
  
  
  }

