
 library(aster2)

 data(test1)
 fred <- asterdata(test1,
     vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
     pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
     code = c(1, 1, 1, 2, 2, 3, 4, 5),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))
 is.validasterdata(fred)

 fred <- asterdata(test1, vars = c("m1", "n1", "n2"), pred = c(0, 1, 1),
     group = c(0, 0, 2), code = c(1, 2, 2),
     families = list("bernoulli", "normal.location.scale"))
 is.validasterdata(fred)

 fred <- asterdata(test1,
     vars = c("m1", "n1", "m2", "b1", "z1", "m3", "p1", "n2"),
     pred = c(0, 1, 0, 3, 4, 0, 6, 1), group = c(0, 0, 1, 0, 0, 3, 0, 2),
     code = c(1, 2, 1, 3, 5, 1, 4, 2),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))
 is.validasterdata(fred)


