library(sprint)
library(RUnit)

# random filled 3 by 4 matrices
mat1 <- matrix(rexp(12), nrow=3)
mat2 <- matrix(rexp(12), nrow=3)

print("about to run pcor")
pcor_result = pcor(mat1, mat2)
cor_result = cor(mat1, mat2)
checkEqualsNumeric(cor_result, pcor_result[,])
print(pcor_result)
