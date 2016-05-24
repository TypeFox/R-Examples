###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)           #
#       Consiglio Nazionale delle Ricerche                   #
#       www.isib.cnr.it                                       #
#                                                             #
#   $Id: zzz.R 388 2015-05-19 19:09:08Z tonig $
#                                                             #
###############################################################

.onAttach <- function(lib, pkg)  {

  packageStartupMessage(sprintf("Loaded dtw v%s. See ?dtw for help, citation(\"dtw\") for use in publication.\n",
            utils::packageDescription("dtw")$Version ) );
      
  ## Register DTW as a distance function into package proxy
  pr_DB$set_entry(FUN=dtwpairdist, names=c("DTW","dtw"), 
                  loop=TRUE, type="metric",
                  description="Dynamic Time Warping",
                  reference="Giorgino T (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: The dtw Package. Journal of Statistical Software, 31 (7), pp. 1--24. <URL: http://www.jstatsoft.org/v31/i07/>.",
                  formula="minimum of sum(x[xw[i]]-y[yw[i]]) over all monotonic xw, yw");

  # invisible()
}
