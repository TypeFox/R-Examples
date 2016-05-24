# Test data from:
#    Batschelet, E (1981). Circular Statistics in Biology.
#    Examples 6.10.1 and 6.10.2, p 126
# 

suppressMessages(library("circular"))
# ?wallraff.test

angles <- circular(c(70, 80, 80, 85, 85, 90, 95, 95, 5, 5, 15, 55, 55, 65, 105, 120, 340), units="degrees", template="geographics")
group <- factor(c(rep("control", 8), rep("experimental", 9)))

homeDir <- 40

# expect:
# W = 2 (in wilcox.test) and p < 0.01 for the dispersion test
# W = 26 (in wilcox.test) and p > 0.05 for the homing test

xn <- angles
wallraff.test(xn, group)

wallraff.test(xn, group, ref=homeDir)
wallraff.test(xn, as.factor(group), ref=homeDir)


xl <- split(xn, group)
wallraff.test(xl, ref=homeDir)
wallraff.test(xl)

xl <- split(xn, group)
names(xl) <- NULL
wallraff.test(xl)

xd <- data.frame(group=group, angles=angles)
wallraff.test(angles ~ group, xd)
