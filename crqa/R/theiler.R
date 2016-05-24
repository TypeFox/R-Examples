## define a theiler window to remove recurrent points along the diagonal of the matrix
## written by Moreno I. Coco (moreno.cocoi@gmail.com)
## S = a sparse recurrent matrix
## tw = the theiler window

# S = rbind(c(0,1,1,1,0), c(0,0,1,1,1), c(1,1,1,0,0),c(1,1,0,0,1),c(1,1,0,1,0))

theiler <- function(S,tw) {

    if (tw > nrow(S)) stop ("Theiler window larger than number of diagonals")

    if (tw > 0){ ## remove the theilers' points
        theilers = vector()
    
        for (t in 1:tw){
            
            if (t == 1){ ## this is the diagonal
                theilers = rbind(theilers, which(row(S) == col(S), arr.ind = TRUE),
                    deparse.level = 0)
            } else {
                
                ## above diagonal
                theilers = rbind(theilers, which(row(S) == col(S) + (t -1), arr.ind = TRUE),
                    deparse.level = 0)
                
                ## below diagonal
                theilers = rbind(theilers, which(row(S) == col(S) - (t -1), arr.ind = TRUE),
                    deparse.level = 0)
                
            }
        }

        S[theilers] = 0
        return (S)
    } else { ## leave the matrix untouched
        return(S)
    }

}

