###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: warpArea.R 311 2013-06-03 17:23:17Z tonig $
#                                                             #
###############################################################


## Compute the (approximate) area between the warping function and the
## diagonal (in unit steps). 


warpArea <- function(d) {

  if(!is.dtw(d))
    stop("dtw object required");
  
  ## rebuild query->templ map, interpolating holes
  ii<-approx(x=d$index1,y=d$index2,1:d$N);

  dg<-seq(from=1,to=d$M,len=d$N);

  ad<-abs(ii$y-dg);
  sum(ad);
     
}




## Exmp:
##   t <- localWarpingStretch(alignment)
##
##   # local warping amount, based on the warping function
##   plot(t)
##
##   # lwa, based on the input position
##   plot(t~alignment$index1)
##
##   # or its pointwise approximation
##   plot(approx(alignment$index1,t,1:alignemnt$N)$y)

##
##   diff(localWarpingStretch) contains +1 for each "insertion" and -1
##    for each "deletion". Since we have max N deletions + M
##    insertions, a reasonable normalization would be to divide
##    sum(abs(diff())) by 2N
##
##   localWarpingCost is the amount of mismatch
##    at each point ("local substitution cost")



## Return how far from the diagonal is each point on the warping
## function. A good normalization factor can be 2 * N, so that maximum
## stretch is 1.

.localWarpingStretch <- function(d) {
    ## The diagonal line
    # dg <- seq(from = 1, to = d$M, len = d$N)

    ## get local copies
    id1 <- d$index1;
    id2 <- d$index2;

    ## remap reference indices to a square alignment
    id2 <- id2*d$N/d$M;

    ## return the local distance from the diagonal
    id1-id2;
}



## Return the local costs along the warping path i.e. d[i[t],j[t]] . A
## reasonable normalization could be to d$distance, so that each
## element would be the fraction of cost accumulated at that step.

.localWarpingCost <- function(d) {
    if(is.null(d$localCostMatrix))
        stop("A dtw object with keep.internals=TRUE is required");

    diag(d$localCostMatrix[d$index1,d$index2]);
}
