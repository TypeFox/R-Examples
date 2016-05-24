"ClassUtil" <-
function (rxy = 0, BR = .5, SR = .5)
 {
 pTP <- BR*SR + rxy*sqrt(BR*(1-BR) * SR*(1-SR))
 pFN <- BR - pTP
 pFP <- SR - pTP
 pTN <- 1 - pTP - pFN - pFP
 sen <- pTP/(pTP+pFN)
 spe <- pTN/(pFP+pTN)
 cd  <- (pTP+pTN)*100
 suc <- pTP/(pTP+pFP)
imp <- (suc - BR)*100

mat <- matrix(rbind(pTP,pFN,pFP,pTN,NA,sen,spe,cd,suc,imp))
 colnames(mat) <- "Probabilities"
 rownames(mat) <- c("True Positives", "False Negatives", "False Positives", 
 "True Negatives","--", "Sensitivity", "Specificity",
 "% of Decisions Correct", "Proportion Selected Succesful",
 "% Improvement over BR")
return(mat)
}

