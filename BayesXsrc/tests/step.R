## BayesX STEP testing
library("BayesXsrc")
step <- run.bayesx("step.prg", verbose = FALSE)
fx1 <- read.table("step_f_x1_pspline.res", header = TRUE)
fx4 <- read.table("step_f_x4_pspline.res", header = TRUE)
print(round(head(fx1), digits = 3))
print(round(head(fx4), digits = 3))
