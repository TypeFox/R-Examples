# $Id: minSigCor.R, v1.2.3 2015/08/05 12:00:00 hsbadr EPS JHU             #
#-------------------------------------------------------------------------#
# This function is a part of HiClimR R package.                           #
#-------------------------------------------------------------------------#
#  HISTORY:                                                               #
#-------------------------------------------------------------------------#
#  Version  |  Date      |  Comment   |  Author          |  Email         #
#-------------------------------------------------------------------------#
#           |  May 1992  |  Original  |  F. Murtagh      |                #
#           |  Dec 1996  |  Modified  |  Ross Ihaka      |                #
#           |  Apr 1998  |  Modified  |  F. Leisch       |                #
#           |  Jun 2000  |  Modified  |  F. Leisch       |                #
#-------------------------------------------------------------------------#
#   1.0.0   |  03/07/14  |  HiClimR   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.1   |  03/08/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.2   |  03/09/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.3   |  03/12/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.4   |  03/14/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.5   |  03/18/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.6   |  03/25/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.0.7   |  03/30/14  |  Hybrid    |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.8   |  05/06/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.0.9   |  05/07/14  |  CRAN      |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.0   |  05/15/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.1   |  07/14/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.2   |  07/26/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.3   |  08/28/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.4   |  09/01/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.5   |  11/12/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.1.6   |  03/01/15  |  GitHub    |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.2.0   |  03/27/15  |  MVC       |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.1   |  05/24/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.2   |  07/21/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.3   |  08/05/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
# COPYRIGHT(C) 2013-2015 Earth and Planetary Sciences (EPS), JHU.         #
#-------------------------------------------------------------------------#
# Function: Minimum significant correlation for a sample size             #
#-------------------------------------------------------------------------#

minSigCor <- function(n = 41, alpha = 0.05, r = seq(0, 1, by = 1e-06)) {

    dof <- n - 2
    Fstat <- r^2 * dof/(1 - r^2)
    p.value <- 1 - pf(Fstat, 1, dof)
    p.value[p.value > alpha] <- NA
    i <- which(p.value == max(p.value, na.rm = TRUE))
    
    RsMin <- list(cor = r[i], p.value = p.value[i])
    
    #gc()
    return(RsMin)
}
