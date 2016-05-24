z <- zoprobit$new()
test <- z$mcunit(minx=0, maxx=2, plot=FALSE)
expect_true(test)