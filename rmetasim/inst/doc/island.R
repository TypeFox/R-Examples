seedMigrationRate <- 1
pollenMigrationRate <- 1
habitats <- 10
stages <- 2

rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=habitats, s=stages, totgen=5000)
rland <- landscape.new.switchparam(rland, mp=0)
rland <- landscape.new.floatparam(rland)

S <- matrix(c(0,0,1,0), nrow = 2, byrow = TRUE)
R <- matrix(c(0,1.1,0,0), nrow = 2, byrow = TRUE)
M <- matrix(c(0,0,0,1), nrow = 2, byrow = TRUE)
rland <- landscape.new.local.demo(rland,S,R,M)

rland <- landscape.new.epoch.island(rland,0,c(0,0),c(0,0),
                          seedMigrationRate, c(1,0), c(1,0),
                          pollenMigrationRate, c(0,1), c(0,1))

rland <- landscape.new.locus(rland,type=2,ploidy=1,transmission=1,numalleles=4,allelesize=100)
for (x in 1:9) {
  rland <- landscape.new.locus(rland,type=2,ploidy=2,transmission=0,numalleles=4,allelesize=100)
}
for (x in 1:10) {
  rland <- landscape.new.locus(rland,type=1,ploidy=2,transmission=0,numalleles=4)
}

# add individuals
rland <- landscape.new.individuals(rland,
                         c(5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0,
                         5, 0))
