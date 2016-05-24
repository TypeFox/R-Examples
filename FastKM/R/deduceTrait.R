deduceTrait <- function(givenDelta, y, trait){

  if( !is(trait,'NULL') ) trait <- tolower(trait)

  if( !givenDelta ) {

    if( !is(trait,'NULL') ) {
      if( trait == 's' ) stop("delta is required for survival traits.",
                              call. = FALSE)
    }

    #--------------------------------------------------------------#
    # If delta is not provided, determine if y is integer.         #
    #--------------------------------------------------------------#
    yBin <- isTRUE(all.equal(round(y,0L),y))

    if( yBin ) {

      #----------------------------------------------------------#
      # If y is integer, determine how many unique values.       #
      #----------------------------------------------------------#
      nyValues <- length(unique(as.integer(round(y,0L))))

      if( nyValues > 2.5 ) {

        #------------------------------------------------------#
        # If > 2 unique values, characterize as continuous.    #
        #------------------------------------------------------#
        if( !is(trait,'NULL') ) {
          if( trait != "c" ) {
            cat("y appears to be integer valued with > 2 levels.\n",
                "Should be continuous. Verify Input.", sep="")
          }
        } else {
          trait <- "Continuous"
          cat("\n", trait, " trait analysis.\n", sep="")
        }

      } else if( (nyValues > 1.5) && (nyValues < 2.5) ) {

        #------------------------------------------------------#
        # If = 2 unique values, characterize as dichotomous.   #
        #------------------------------------------------------#
        if( !is(trait, 'NULL') ) {
          if( trait != "d" ) {
            cat("y appears to be integer valued with 2 levels.\n",
                "Should be dichotomous. Verify Input.", sep="")
          }
        } else {
          trait <- "Dichotomous"
          cat("\n", trait, " trait analysis.\n", sep="")
        }

      } else {

        #------------------------------------------------------#
        # If < 2 unique values, stop with error.               #
        #------------------------------------------------------#
        stop("Unable to identify trait type.", call. = FALSE)

      }
    } else {

      #----------------------------------------------------------#
      # If y is non-integer, characterize as continuous.         #
      #----------------------------------------------------------#
      if( !is(trait,'NULL') ) {
        if( trait != "c" ) {
          cat("y appears to be continuous. Verify input.\n")
        }
      } else {
        trait <- "Continuous"
        cat("\n", trait, " trait analysis.\n", sep="")
      }

    }
  } else {

    #--------------------------------------------------------------#
    # If delta is provided, characterize trait as survival.        #
    #--------------------------------------------------------------#
    if( !is(trait,'NULL') ) {
      if( trait != "s" ) {
        stop("delta is provided. y must be of type survival.\n")
      }
    } else {
      trait <- "Survival"
      cat("\n", trait, " trait analysis.\n", sep="")
    }

  }

  trait <- substr(tolower(trait),1,1)

  return(trait)


}
