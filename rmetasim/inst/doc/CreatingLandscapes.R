## ------------------------------------------------------------------------
library(rmetasim)
rland <- landscape.new.empty() #creates a skeleton of a landscape
rland <- landscape.new.intparam(rland, h=2, s=2) 

names(rland$intparam)

## ------------------------------------------------------------------------
rland <- landscape.new.switchparam(rland)
names(rland$switchparam)

## ------------------------------------------------------------------------
rland <- landscape.new.floatparam(rland)
names(rland$floatparam)

## ------------------------------------------------------------------------

S <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
R <- matrix(c(0, 1.1, 0, 0), nrow = 2)
M <- matrix(c(0, 0, 0, 1), nrow = 2)
rland <- landscape.new.local.demo(rland,S,R,M)

S <- matrix(c(rep(0,4),
              rep(0,4),
              rep(0,4),
              rep(0,4)), nrow = 4)

R <- matrix(c(0,0,0,1,
              0,0,0,0,
              0,1,0,0,
              0,0,0,0), byrow=T, nrow = 4)

M <- matrix(c(0,0,0,0,
              0,0,0,.1,
              0,0,0,0,
              0,.1,0,0), nrow = 4)

rland <- landscape.new.epoch(rland,S=S,R=R,M=M,extinct=c(.01,0),carry=c(100,200))

names(rland$demography)
names(rland$demography$localdem[[1]])
names(rland$demography$localdemK[[1]])
names(rland$demography$epochs[[1]])

## ------------------------------------------------------------------------
rland <- landscape.new.locus(rland,type=0,ploidy=2,mutationrate=0.001,
                   transmission=0,numalleles=5)
rland <- landscape.new.locus(rland,type=1,ploidy=1,mutationrate=0.005,
                   numalleles=3,frequencies=c(.2,.2,.6))
rland <- landscape.new.locus(rland,type=2,ploidy=2,mutationrate=0.007,
                   transmission=0,numalleles=6,allelesize=75)


names(rland$loci[[1]])
names(rland$loci[[1]]$alleles[[1]])

## ------------------------------------------------------------------------
rland <- landscape.new.individuals(rland,c(10,15,20,8))
print(rland$individuals[1,])

## ------------------------------------------------------------------------
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h = 2, s = 2)
rland <- landscape.new.switchparam(rland, mp = 0)
rland <- landscape.new.floatparam(rland)

S <- matrix(c(0, 0, 1, 0), byrow = TRUE, nrow = 2)
R <- matrix(c(0, 1.1, 0, 0), byrow = TRUE, nrow = 2)
M <- matrix(c(0, 0, 0, 1), byrow = TRUE, nrow = 2)

rland <- landscape.new.local.demo(rland, S, R, M)

S <- matrix(rep(0, 16), nrow = 4)
R <- matrix(rep(0, 16), nrow = 4)
M <- matrix(rep(0, 16), nrow = 4)

rland <- landscape.new.epoch(rland, S = S, R = R, M = M, 
                             carry = c(1000, 1000))

rland <- landscape.new.locus(rland, type = 0, ploidy = 2, 
                             mutationrate = 0.001, transmission = 0, numalleles = 5)
rland <- landscape.new.locus(rland, type = 1, ploidy = 1, 
                             mutationrate = 0.005, numalleles = 3, frequencies = c(0.2, 
                                                                       0.2, 0.6))
rland <- landscape.new.locus(rland, type = 2, ploidy = 2, 
                             mutationrate = 0.007, transmission = 0, numalleles = 6, 
                             allelesize = 75)
rland <- landscape.new.individuals(rland, c(50, 0, 50, 0))

## ------------------------------------------------------------------------
library(magrittr)
#first set up the matrices for local demographies
S <- matrix(c(0, 0, 1, 0), byrow = TRUE, nrow = 2)
R <- matrix(c(0, 1.1, 0, 0), byrow = TRUE, nrow = 2)
M <- matrix(c(0, 0, 0, 1), byrow = TRUE, nrow = 2)

#and epochs 
S.epoch <- matrix(rep(0, 16), nrow = 4)
R.epoch <- matrix(rep(0, 16), nrow = 4)
M.epoch <- matrix(rep(0, 16), nrow = 4)

##now create the landscape
rland <- landscape.new.empty() %>% 
    landscape.new.intparam(h=2,s=2) %>% 
    landscape.new.switchparam(mp=0) %>%
    landscape.new.floatparam() %>%
    landscape.new.local.demo( S, R, M) %>%
    landscape.new.epoch(S = S.epoch, R = R.epoch, M = M.epoch, 
                        carry = c(1000, 1000)) %>%
    landscape.new.locus(type = 0, ploidy = 2, 
                        mutationrate = 0.001, transmission = 0, numalleles = 5) %>%
    landscape.new.locus(type = 1, ploidy = 1, 
                        mutationrate = 0.005, numalleles = 3, frequencies = c(0.2, 
                                                                  0.2, 0.6)) %>%
    landscape.new.locus(type = 2, ploidy = 2, 
                        mutationrate = 0.007, transmission = 0, numalleles = 6, 
                        allelesize = 75) %>%
    landscape.new.individuals(c(50, 0, 50, 0))
           

## ------------------------------------------------------------------------
rland <- landscape.simulate(rland,10)
landscape.amova(rland)

## ------------------------------------------------------------------------
rland %>% landscape.simulate(10) %>% landscape.amova()

