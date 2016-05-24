

tspot.criteria.10spot <- function(tspot.obj){

    
    # Set up the results vector
    result <- rep(NA, times = length(tspot.obj$nil)) 

    # Identify the maximum of Panel A - Nil and Panel B - Nil
    panel.max <- pmax((tspot.obj$panel.a - tspot.obj$nil),
                      (tspot.obj$panel.b - tspot.obj$nil)
    )

    # For ease of reading criteria, set up a flag for each of the TSPOT
    # criteria
    # Nil is high
    high.nil <- tspot.obj$nil > 10

    # Antigen is in borderline territory
    border.spots <- panel.max %in% 5:9

    # Antigen is in positive territory
    pos.spots <- panel.max >= 10

    # Mitogen is low
    low.mito <- tspot.obj$mito < 20


    # Positive
    result[border.spots %in% FALSE &
           pos.spots %in% TRUE &
           high.nil %in% FALSE &
           low.mito %in% c(TRUE, FALSE)] <- "Positive"

    # Borderline
    result[border.spots %in% TRUE &
           pos.spots %in% FALSE &
           high.nil %in% FALSE &
           low.mito %in% c(TRUE, FALSE)] <- "Borderline"


    # Negative
    result[border.spots %in% FALSE &
           pos.spots %in% FALSE &
           high.nil %in% FALSE &
           low.mito %in% FALSE] <- "Negative"


    # High nil
    result[border.spots %in% c(TRUE, FALSE) &
           pos.spots %in% c(TRUE, FALSE) &
           high.nil %in% TRUE &
           low.mito %in% c(TRUE, FALSE)] <- "Invalid - high nil"


    # Low mito
    result[border.spots %in% FALSE &
           pos.spots %in% FALSE &
           high.nil %in% FALSE &
           low.mito %in% TRUE] <- "Invalid - low mitogen"




    return(result)
}

