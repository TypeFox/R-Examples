
cell(3)
cell(2,3)

## define possible values of 2 categorical design variable
x.space <- cell(1,2)
x.space[1,1] <- list(seq(10,100,10))
x.space[1,2] <- list(seq(10,300,10))
x.space
x.space[1,1]
x.space[1,2]
