
 library(aster2)

 data(test1)

 fred <- asterdata(test1,
     vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
     pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
     code = c(1, 1, 1, 2, 2, 3, 4, 5),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))

 theta.star <- aster2:::starting(fred)

 idx.minus.one <- as.character(fred$redata$varb) == "n2"

 all(theta.star[idx.minus.one] == (- 1))
 all(theta.star[! idx.minus.one] == 0)

