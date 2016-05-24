se###
###
### Ramsey & Silverman (2006) Functional Data Analysis, 2nd ed. (Springer)
###
### ch. 17.  Derivatives and functional linear models 
###
library(fda)
##
## Section 17.1.  Introduction
##
# Derivatives can be as important as a function differentiated

##
## Section 17.2.  Oil refinery data 
##
#. p. 298, Figure 17.1.  Reflux flow and Tray 47 level in a refinery 

str(refinery)
sapply(refinery, class)

with(growth, matplot(age, hgtf[, sel], type="b"))

refOrder <- 4
(refKnots <- c(0, 67/2, 67, 67, 67, 67+(1:4)*(193-67)/4))
(ref.nKnots <- length(refKnots))
# p. 299:  "These knot choices imply a total of eleven basis functions."  
(ref.nBasis <- ref.nKnots+refOrder-2)

# I should work fda-ch14 before continuing here.  
                   
