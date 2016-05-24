## roundConstantSum
## 
## Ian Kopacka
## 2011-07-05
##
## Inputparameter:
##    numVec...Numerischer Vektor. Die Summe der Elemente sollte eine
##                   ganze Zahl sein. Falls nicht, wird mit einem gerundeten Wert
##                   gerechnet und wahlweise eine Warnung augegeben.
##    output.........0 oder 1. Bei 0 wird keine Warnung ausgegeben falls 
##                   'sum(numVec)' keine ganze Zahl ist. Bei 1 (=default)
##                   wird eine Warnung auf dem Bildschirm ausgegeben.
## Outputparameter: numerischer Vektor der selben Laenge wie 'numVec'. Die 
## Elemente sind ganze Zahlen, und 
##           sum(roundConstantSum(numVec)) == sum(numVec)

roundConstantSum <- function(numVec, output = 1){
    Summe <- sum(numVec)
    if(Summe - round(Summe) > 0){
        if(output == 1){
            warnung <- paste("WARNUNG - Sum of 'numVec' is not an integer: ", Summe, 
                "\n Proceeding with rounded value: ", round(Summe), "\n", sep = "")
            cat(warnung)
        }
        Summe <- round(Summe)    
    }
    gerundet <- round(numVec)
    ## Es wurde zu oft abgerundet:
    if(sum(gerundet) < Summe){
        ## Wie oft wurde faelschlicherweise abgerundet:
        anzFehler <- Summe - sum(gerundet)
        ## Konzept: Suche von den Zahlen, die abgerundet wurden 
        ## die mit den groessten Nachkommastellen heraus und 
        ## runde sie auf:
        rest <- numVec - gerundet      
        names(rest) <- seq(along = rest)
        rest <- rest[order(numVec, decreasing = TRUE)]
        rest <- sort(rest, decreasing = TRUE)
        index <- as.numeric(names(rest)[1:anzFehler])
        gerundet[index] <- gerundet[index] + 1   
        return(gerundet)
    }
    ## Es wurde zu oft aufgerundet:
    if(sum(gerundet) > Summe){
        ## Wie oft wurde faelschlicherweise aufgerundet:
        anzFehler <- sum(gerundet) - Summe
        ## Konzept: Suche von den Zahlen die aufgerundet wurden
        ## die mit den kleinsten Nachkommastellen heraus und
        ## runde sie ab:
        rest <- numVec - gerundet
        names(rest) <- seq(along = rest)
        rest <- rest[order(numVec, decreasing = FALSE)]
        rest <- sort(rest, decreasing = FALSE)
        index <- as.numeric(names(rest)[1:anzFehler])
        gerundet[index] <- gerundet[index] - 1   
        return(gerundet)
    }    
    return(gerundet)
}


