
library(random)

X <- randomNumbers()
dim(X)                                  # cannot easily summary stats or content
#min(X)
#max(X)

X <- randomSequence()                   # here we can test min and max as it is just a shuffle
dim(X)
min(X)
max(X)

X <- randomStrings()
dim(X)
