## ------------------------------------------------------------------------
library(rmetasim)
rland <- landscape.new.example()
rland <- landscape.simulate(rland,10)
landscape.amova(rland)

## ------------------------------------------------------------------------
rland <- landscape.new.example()
for (gen in 1:10)
    {
        rland <- landscape.extinct(rland)
        rland <- landscape.reproduce(rland)
        rland <- landscape.survive(rland)
        rland <- landscape.carry(rland)
        rland <- landscape.advance(rland)
    }
landscape.amova(rland)

## ------------------------------------------------------------------------
library(magrittr)
rland <- landscape.new.example()
for (gen in 1:10)
    {
        rland <- rland %>% landscape.extinct() %>%  landscape.reproduce() %>% 
            landscape.survive() %>% landscape.carry()%>%landscape.advance()
    }
landscape.amova(rland)

## ------------------------------------------------------------------------
create.land <- function()
    {
### CREATE A LANDSCAPE
###first set up the matrices for local demographies
        S <- matrix(c(0, 0, 1, 0), byrow = TRUE, nrow = 2)
        R <- matrix(c(0, 1.1, 0, 0), byrow = TRUE, nrow = 2)
        M <- matrix(c(0, 0, 0, 1), byrow = TRUE, nrow = 2)
        
                                        #and epochs 
        S.epoch <- matrix(rep(0, 36), nrow = 6)
        R.epoch <- matrix(rep(0, 36), nrow = 6)
        M.epoch <- matrix(rep(0, 36), nrow = 6)
        
        
        
        ##now create the landscape
        landscape.new.empty() %>% 
            landscape.new.intparam(h=3,s=2) %>% 
                landscape.new.switchparam(mp=0) %>%
                    landscape.new.floatparam() %>%
                        landscape.new.local.demo( S, R, M) %>%
                            landscape.new.epoch(S = S.epoch, R = R.epoch, M = M.epoch, 
                              carry = c(1000, 1000, 1000)) %>%
                                landscape.new.locus(type = 1, ploidy = 2, 
                                  mutationrate = 0.001, transmission = 0, numalleles = 5) %>%
                                landscape.new.locus(type = 0, ploidy = 1, 
                                  mutationrate = 0.005, numalleles = 3, frequencies = c(0.2, 0.2, 0.6)) %>%
                                landscape.new.locus(type = 2, ploidy = 1, 
                                  mutationrate = 0.007, transmission = 1, numalleles = 3, 
                                  allelesize = 75) %>%
                                landscape.new.individuals(rep(c(500,500),3))
    }

## ------------------------------------------------------------------------
rland <- create.land()
retval <- matrix(0,ncol=4,nrow=21) #store the results col1 = gen, cols2-3 PhiST
retval[1,] <- c(0,landscape.amova(rland)) #before any evolution occurs. Both populations from same source
for (i in 2:21)
    {
        rland <- landscape.simulate(rland,5)
        retval[i,] <- c(rland$intparam$currentgen,landscape.amova(rland)) 
    }
###the rest is just plotting the matrix of PhiSTs
###you can use any kind of graphics, of course.  using ggplot2 here.
colnames(retval) <- c("gen","Loc1","Loc2","Loc3")
library(ggplot2)
library(reshape2)
pl <- ggplot(melt(as.data.frame(retval),measure.vars=c("Loc1","Loc2","Loc3")),
             aes(x=gen,y=value,group=variable,color=variable))
pl <- pl + labs(y="PhiST",x="Year",color="Locus")
pl <- pl + geom_point() + geom_smooth()
pl

## ------------------------------------------------------------------------
number.reps=3
retlst <- lapply(1:number.reps,function(x) 
                 {
                     rland <- create.land()
                     retval <- matrix(0,ncol=5,nrow=21) 
                     retval[1,] <- c(0,landscape.amova(rland),x)
                     for (i in 2:21)
                         {
                             rland <- landscape.simulate(rland,5)
                             retval[i,] <- c(rland$intparam$currentgen,landscape.amova(rland),x) 
                         }
                     colnames(retval) <- c("gen","Loc1","Loc2","Loc3","replicate")
                     retval
                 })
retdf <- as.data.frame(do.call(rbind,retlst))

pl <- ggplot(melt(retdf,measure.vars=c("Loc1","Loc2","Loc3")),
             aes(x=gen,y=value,group=variable,color=variable))
pl <- pl + labs(y="PhiST",x="Year",color="Locus")
pl <- pl + geom_point() + geom_smooth()
pl


