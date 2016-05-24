#
# these are functions that update elements of a landscape (often requested)
# AES 1/1/05
#
landscape.modify.epoch <- function(rland,epoch=1,S=NULL,R=NULL,M=NULL,epochprob=NULL,startgen=NULL,extinct=NULL,carry=NULL,localprob=NULL)
  {
    if (!is.null(rland$demography$epochs[[epoch]]))
        {
          if (!is.null(S))
            {
              rland$demography$epochs[[epoch]]$S <- S
            }
          if (!is.null(R))
            {
              rland$demography$epochs[[epoch]]$R <- R
            }
          if (!is.null(M))
            {
              rland$demography$epochs[[epoch]]$M <- M
            }
          if (!is.null(epochprob))
            {
              rland$demography$epochs[[epoch]]$RndChooseProb <- epochprob
            }
          if (!is.null(startgen))
            {
              rland$demography$epochs[[epoch]]$Startgen <- startgen
            }
          if (!is.null(extinct))
            {
              rland$demography$epochs[[epoch]]$Extinct <- extinct
            }
          if (!is.null(carry))
            {
              rland$demography$epochs[[epoch]]$Carry <- carry
            }
          if (!is.null(localprob))
            {
              rland$demography$epochs[[epoch]]$Localprob <- localprob
            }
        }
    rland
  }
