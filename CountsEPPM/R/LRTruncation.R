LRTruncation <-
function (probability,ltvalue,utvalue) {
    nmax1 <- length(probability)
    rev.probability <- rep(NA,nmax1)
    t.prob <- sum(probability)
# rounding error can cause the sum of the probabilities to <0 or >1
# so rounded to 10 decimal places                                        
    round.t.prob <- round(t.prob,digits=10)
    if ((is.finite(t.prob))==TRUE) { 
       if ((round.t.prob>0) & (round.t.prob<=1)) { 
          wks <- sum((probability==0))
          wks1 <- 1
          wks2 <- nmax1
          if (wks<nmax1) { x <- 0 
             if (is.na(utvalue)==FALSE) {    x <- 1 - t.prob
                                          wks2 <- utvalue
                if (nmax1>utvalue) { x <- x + sum(probability[utvalue:nmax1]) } }
             if (is.na(ltvalue)==FALSE) { wks1 <- ltvalue + 1
                                             x <- x + sum(probability[1:wks1])
                                          wks1 <- wks1 + 1 }
             if ((x>0) & (x<1)) { rev.probability[wks1:wks2] <- 
                                               probability[wks1:wks2] / (1 - x)
                                } else { if (x<=0) { 
                                            rev.probability <- probability } }
                                      } } # end of ((t.prob>0) & (t.prob<1)) 
                                   } # end of is.finite(t.prob)
return(rev.probability) }
