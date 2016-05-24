library(FrF2)
### test programs for FrF2
### pb functionality
### combination of everything with everything is not included, 
###    but presumably also not needed
pb(8,randomize=FALSE)
pb(8,nfactors=4,randomize=FALSE)
## n12 with and without Taguchi
pb(12,n12.taguchi=FALSE,randomize=FALSE, default.levels=c("-","+"))
   ## prevent error that both default.levels are identical ?
   try(pb(12,n12.taguchi=FALSE,randomize=FALSE, default.levels=c("-","-")))
pb(12,n12.taguchi=FALSE,randomize=FALSE, replications=2, repeat.only=TRUE)
pb(12,n12.taguchi=FALSE,randomize=FALSE, replications=2, repeat.only=FALSE)
set.seed(98776)
pb(12,n12.taguchi=FALSE,replications=2, repeat.only=FALSE)
pb(12,n12.taguchi=FALSE,randomize=FALSE, factor.names=list(eins=c(1,2),zwei=c(43,87),drei=c("alt","neu"),vier="",fuenf="",sechs=""))
pb(12,n12.taguchi=TRUE,randomize=FALSE)
## box-tyssedal
pb(16,nfactors=12,randomize=FALSE)
pb(16,nfactors=13,randomize=FALSE)
pb(16,nfactors=15,randomize=FALSE)
## doubling
pb(32, randomize=FALSE)
## special cases
pb(28, randomize=FALSE)
pb(52, randomize=FALSE)
pb(76, randomize=FALSE)
pb(92, randomize=FALSE)
