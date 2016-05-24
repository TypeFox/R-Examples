# Test data from:
#    Circular statistics in biology, Batschelet, E (1981)
#    §6.2, p99
# 

suppressMessages(library("circular"))
# ?watson.williams.test

angles <- circular( c(rep(c(-20, -10, 0), c(1,7,2)), rep(c(-10, 0, 10, 20), c(3,3,3,1))), units="degrees", template="geographics")
group <- factor(rep(c("exp", "control"), each=10))

# expect this:
# F = 8.7329, df1 = 1, df2 = 18, p-value = 0.003108
# mean of control     mean of exp 
#        1.988969       -9.000615 

# Test interfaces
xn <- angles
watson.williams.test(xn, group)

xl <- split(xn, group)
watson.williams.test(xl)

xl <- split(xn, group)
names(xl) <- NULL
watson.williams.test(xl)

xd <- data.frame(group=group, angles=angles)
watson.williams.test(angles ~ group, xd)

# Test the influence of ordering the groups
id <- c(9, 8, 7, 4, 6, 5, 12, 18, 10, 17, 1, 19, 3, 20, 2, 16, 15, 14, 13, 11)
angles <- angles[id]
group <- group[id]

xn <- angles
watson.williams.test(xn, group)
xl <- split(xn, group)
watson.williams.test(xl)
xd <- data.frame(group=group, angles=angles)
watson.williams.test(angles ~ group, xd)

# Test NAs
angles[length(angles)+1] <- NA
levels(group) <- c("exp", "control", "bar")
group[length(group)+1] <- "bar"
xn <- angles
watson.williams.test(xn, group)
