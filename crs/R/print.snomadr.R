#
# File:   print.snomadr.R
# Author: Zhenghua Nie
# Date:   Mon 16 May 2011
#
# We use ipoptr developed by Jelmer Ypma as the prototype of this package.
# Some code is copied and edited from ipoptr. 
# Please reference the license of ipoptr.
#
# This function prints some basic output of a snomadr 
# object. The information is only available after it 
# has been solved.
#
# Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
# This code is published under GNU GENERAL PUBLIC LICENSE.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,  or
# (at your option) any later version.
#      
# This program is distributed WITHOUT ANY WARRANTY. See the
# GNU General Public License for more details.
#           
# If you do not have a copy of the GNU General Public License,  
# write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

print.snomadr <- function(x, show.controls=TRUE, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "", fill=TRUE)

  if (!is.null(x$call$information)) {return(NULL);}
  else
  {
      cat( unlist(strsplit(paste( "nomad solver status:", x$status, "(", x$message, ")\n" ),' ')), fill=TRUE )
      cat( paste( "Number of blackbox evaluations.....:", x$bbe, "\n" ) )
      cat( paste( "Number of iterations...............:", x$iterations, "\n" ) )

      # if show.controls is TRUE or FALSE, show all or none of the controls
      if ( is.logical( show.controls ) ) {
          # show all control variables
          if ( show.controls ) {
              controls.indices = 1:length(x$solution)
          }
      }

      # if show.controls is a vector with indices, rename this vector
      # and define show.controls as TRUE
      if ( is.numeric( show.controls ) ) {
          controls.indices = show.controls
          show.controls = TRUE
      }

      # if solved successfully
      if ( x$status == 8 || x$status == 9 || x$status == 13 ) {
          cat( paste( "Optimal value of objective function: ", x$objective, "\n" ) )
          if ( show.controls ) {
              cat( "Optimal value of controls..........: " )
              cat( x$solution[ controls.indices ], fill=TRUE)
              cat("\n")
          }
      } else {
          cat( paste( "Current value of objective function: ", x$objective, "\n" ) )
          if ( show.controls ) {
              cat( "Current value of controls..........: " )
              cat( x$solution[ controls.indices ], fill=TRUE )
              cat("\n")
          }
      }
      cat("\n")
  }
}
