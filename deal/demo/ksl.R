library(deal)  ## invoke DEAL

data(ksl)      ## read data (included in DEAL)

## specify prior network
ksl.nw  <- network(ksl)

## make joint prior distribution
ksl.prior <- jointprior(ksl.nw,64)

## ban arrows towards Sex and Year
mybanlist <- matrix(c(5,5,6,6,7,7,9,
                    8,9,8,9,8,9,8),ncol=2)
banlist(ksl.nw) <- mybanlist                  

## learn the initial network
ksl.nw <- getnetwork(learn(ksl.nw,ksl,ksl.prior))

## Do structural search 
ksl.search <- autosearch(ksl.nw,ksl,ksl.prior,
                      trace=TRUE)
 
## perturb 'thebest' and rerun search twice.
ksl.heuristic <- heuristic(getnetwork(ksl.search),ksl,
                         ksl.prior,
                         restart=2,degree=10,
                         trace=TRUE,trylist=gettrylist(ksl.search))

thebest2 <- getnetwork(ksl.heuristic)

## Run this to transfer network to .net format
## (read by e.g. Hugin, www.hugin.com)
# savenet(thebest2, file("ksl.net"))
