library(MCDA)

style <- t(matrix(c(1,0.25,4,1/6,4,1,4,0.25,0.25,0.25,1,0.2,6,4,5,1),
                  nrow=4,ncol=4))

colnames(style) = c("Corsa","Clio","Fiesta","Sandero")
rownames(style) = c("Corsa","Clio","Fiesta","Sandero")

reliability <- t(matrix(c(1,2,5,1,0.5,1,3,2,0.2,1/3,1,0.25,1,0.5,4,1),
                        nrow=4,ncol=4))

colnames(reliability) = c("Corsa","Clio","Fiesta","Sandero")
rownames(reliability) = c("Corsa","Clio","Fiesta","Sandero")

fuel <- t(matrix(c(1,2,4,1,0.5,1,3,2,0.25,1/3,1,0.2,1,0.5,5,1),nrow=4,ncol=4))

colnames(fuel) = c("Corsa","Clio","Fiesta","Sandero")
rownames(fuel) = c("Corsa","Clio","Fiesta","Sandero")

alternativesPairwiseComparisonsList <- list(style=style, 
                                            reliability=reliability, 
                                            fuel=fuel)

criteriaWeightsPairwiseComparisons <- t(matrix(c(1,0.5,3,2,1,4,1/3,0.25,1),
                                               nrow=3,ncol=3))
colnames(criteriaWeightsPairwiseComparisons) = c("style","reliability","fuel")
rownames(criteriaWeightsPairwiseComparisons) = c("style","reliability","fuel")

overall1 <- AHP(criteriaWeightsPairwiseComparisons, 
                alternativesPairwiseComparisonsList)

s<-structure(c(0.2922, 0.2767, 0.0697, 
               0.3613), .Names = c("Corsa", "Clio", "Fiesta", "Sandero"
               ))

stopifnot(round(overall1,4) == s)