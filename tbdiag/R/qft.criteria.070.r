


################################################################################
# The default Cellestis USA criteria with an increased TB-Nil threshold

qft.criteria.070 <- function(qft.obj){
    # This method calculates the QFT interpretation based on
    # Cellestis' American criteria with a 0.70 IU/mL cutoff instead of the
    # default 0.35 IU/mL cutoff

    # Floating point comparisons can be a problem here.
    # Instead of >=, define a small value and add it to the number 
    # being compared.
    # In essence, convert any left-hand value that's truly equal
    # to the right-hand # value into one that is *greater* 
    # than the right-hand value.
    # This is the tolerance value used by all.equal().
    tol <- .Machine$double.eps ^ 0.5

    
################################################################################
    # Set up the results vector
    result <- rep(NA, times = length(qft.obj$nil)) 



################################################################################
    # For ease of reading criteria, set up a flag for each of the QFT
    # flowchart criteria
    # TB - Nil >= 0.70
    tbnildiff.abs <- (qft.obj$tb - qft.obj$nil) + tol >= 0.70

    # TB - Nil >= 0.25*Nil
    tbnildiff.rel <- (qft.obj$tb - qft.obj$nil) + tol >= (.25 * qft.obj$nil)

    # Nil > 8
    highnil <- qft.obj$nil > 8.0

    # Mito - Nil < 0.50
    mitonildiff <- (qft.obj$mito - qft.obj$nil) + tol < 0.5



################################################################################
    # Compute the results
    
    # Positive
    result[tbnildiff.abs %in% TRUE &
           tbnildiff.rel %in% TRUE &
           highnil %in% FALSE &
           mitonildiff %in% c(TRUE, FALSE)] <- "Positive"

    # Negative
    # TB - Nil is less than 0.70 OR TB - Nil is greater than 0.70, but less
    # than 0.25% of Nil
    result[(tbnildiff.abs %in% FALSE |
           (tbnildiff.abs %in% TRUE & tbnildiff.rel %in% FALSE)) &
           highnil %in% FALSE &
           mitonildiff %in% FALSE] <- "Negative"

    # Indeterminate - high nil
    result[tbnildiff.abs %in% c(TRUE, FALSE) &
           tbnildiff.rel %in% c(TRUE, FALSE) &
           highnil %in% TRUE &
           mitonildiff %in% c(TRUE, FALSE)] <- "Indeterminate - high nil"

    # Indeterminate - mito too close to nil
    result[(tbnildiff.abs %in% FALSE |
           (tbnildiff.abs %in% TRUE & tbnildiff.rel %in% FALSE)) &
           highnil %in% FALSE &
           mitonildiff %in% TRUE] <- "Indeterminate - mitogen too close to nil"



    return(result)

}

