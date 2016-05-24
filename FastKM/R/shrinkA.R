shrinkA <- function(A, aEV, bEV, nZa, nZb, nX) {

  if( nZa == 1L && nZb == 1L ) {
    #--------------------------------------------------------------#
    # If both PVs are of length 1, warn user that a stable         #
    # solution could not be found; exit.                           #
    #--------------------------------------------------------------#
    warning("Unable to identify a stable augmentation matrix.", 
            call.=FALSE)
    return(NULL)

  } else if( nZa == 1L ) {
    #--------------------------------------------------------------#
    # If the PV of Za is of length 1, remove a column from Zb.     #
    #--------------------------------------------------------------#
    A <- A[,-{nX + nZa + nZb}]
    nZb <- nZb - 1L

  } else if( nZb == 1L ) {
    #--------------------------------------------------------------#
    # If the PV of Zb is of length 1, remove a column from Za.     #
    #--------------------------------------------------------------#
    A <- A[,-{nX+nZa}]
    nZa <- nZa - 1L

  } else {
    if( isTRUE(all.equal(aEV[nZa-1L], bEV[nZb-1L])) ) {
      #----------------------------------------------------------#
      # If the next-to-last PV of Za and the next-to-last PV     #
      # of Zb are the same, determine which PV has the most      #
      # elements remaining.                                      #
      #----------------------------------------------------------#
      if( nZa == nZb ) {
        #------------------------------------------------------#
        # If aEV and bEV are the same length, remove last      #
        # element from Za.                                     #
        #------------------------------------------------------#
        A <- A[,-{nX + nZa}]
        nZa <- nZa - 1L

      } else if( nZa > nZb ) {
        #------------------------------------------------------#
        # If aEV longest, remove last element from Za.         #
        #------------------------------------------------------#
        A <- A[,-{nX + nZa}]
        nZa <- nZa - 1L

      } else if( nZa < nZb ) {
        #------------------------------------------------------#
        # If bEV longest remove last element from Zb.          #
        #------------------------------------------------------#
        A <- A[,-{nX + nZa + nZb}]
        nZb <- nZb - 1L

      }
    } else if( aEV[nZa-1L] < bEV[nZb-1L] ) {
      #----------------------------------------------------------#
      # If next-to-last tPV of Zb is largest remove last         #
      # element from Zb.                                         #
      #----------------------------------------------------------#
      A <- A[,-{nX + nZa + nZb}]
      nZb <- nZb - 1L

    } else if( aEV[nZa-1L] > bEV[nZb-1L] ) {
      #----------------------------------------------------------#
      # If next-to-last PV of Za is largest remove an element    #
      # from Za.                                                 #
      #----------------------------------------------------------#
      A <- A[,-{nX + nZa}]
      nZa <- nZa - 1L
    }
  }

  return(list(  "A" = A, 
              "nZa" = nZa, 
              "nZb" = nZb))
}

