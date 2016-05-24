###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: dtw-internal.R 267 2012-08-12 14:37:26Z tonig $
#                                                             #
###############################################################

##
## $Id: dtw-internal.R 267 2012-08-12 14:37:26Z tonig $
##

## Internal functions for the dtw package.
## Not to be used by the user.


## Function applying dtw and only returning the
## distance
dtwpairdist <- function(...) {
  dtw(distance.only=TRUE,...)$distance;
}
