####################################
# Software Testing for mpm package #
####################################

# TODO add graphicsqc-based tests

library(mpm)

# Weighted spectral map analysis
data(Golub) # Gene expression data of leukemia patients
data(Golub.grp) # Pathological classes coded as 1, 2, 3
r.sma <- mpm(Golub[,1:39], row.weight = "mean", col.weight = "mean")

# simple plot
plot(r.sma, label.tol = 20, scale = "uvc",
    col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5)

# change default of main argument
plot(r.sma, label.tol = 20, scale = "uvc",
    col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5, 
    main = "test main")

# change default of sub argument 
plot(r.sma, label.tol = 20, scale = "uvc",
    col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5, sub = "")

# omit sample names and add smooth scatter 
plot(r.sma, label.tol = 20, scale = "uvc",
    col.group = (Golub.grp)[1:38], zoom = c(1,1.2), col.size = 5, sub = "",
    sampleNames = FALSE, do.smoothScatter = TRUE)

# summary and export
export(summary(r.sma))
