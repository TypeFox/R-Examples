# test.earth.R
# The intention here is to check for porting problems by
# building a simple model.
# For more comprehensive tests see earth\inst\slowtests.
library(earth)
options(digits=3)
data(trees)
earth.mod <- earth(Volume~., data=trees)
print(summary(earth.mod))
