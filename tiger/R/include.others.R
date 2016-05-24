include.others <- function(selected, center, stats, best=FALSE){
    above <- sign(stats[3,selected] - center) == 1
    if(best){
        if(above){
           other <- stats[3,] > center &
                    stats[3,] <   #median
                    stats[4,selected] #upper hinge of selected 
        } else { #below
           other <- stats[3,] < center &
                    stats[3,] >   #median
                    stats[2,selected] #lower hinge of selected 
        }
    } else {
        if(above){
           other <- stats[3,] >   #median
                    stats[2,selected] & #lower hinge of selected 
                    stats[2,] > center  #value for no-error does not
                                        #fall within interquartile range
        } else { #below
           other <- stats[3,] <   #median
                    stats[4,selected] & #upper hinge of selected 
                    stats[4,] < center  #value for no-error does not
                                        #fall within interquartile range

        }
    }
    other <- which(other)
    return(other)
}
