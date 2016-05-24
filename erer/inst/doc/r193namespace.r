# A. packages avaiable
options(width = 73); p1 <- available.packages(); str(p1)
p1[c("apt", "erer", "xlsx"), c(1, 2, 4)]

# B. packages installed on a computer
p2 <- installed.packages(); str(p2)
p2[c("apt", "erer", "xlsx"), c(1, 3, 5)]

# C. Packages loaded and attached
# C1. At the beginning of an R session
sa <- sessionInfo(); sa
listn  # error: object 'listn' not found

# C2. Loading "apt" library only
loadNamespace(package = "apt")
sb <- sessionInfo(); sb
listn        # error: object 'listn' not found
erer::listn  # work well

# C3. Attaching "apt" library
library(apt)
sc <- sessionInfo(); sc
listn        # work well
erer::listn  # work well