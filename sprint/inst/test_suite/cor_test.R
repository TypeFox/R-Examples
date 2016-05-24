library("sprint")

my.matrix <- matrix(rnorm(10000,9,1.7), nrow=70000, ncol=25)
genecor <- cor( t(my.matrix) )
print(genecor)

