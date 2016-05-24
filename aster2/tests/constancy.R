
 library(aster2)

 ##### parm.type = "theta"

 data(echinacea)
 cmat <- constancy(echinacea, parm.type = "theta")
 nrow(cmat) == 0

 data(test1)
 fred <- asterdata(test1,
     vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
     pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
     code = c(1, 1, 1, 2, 2, 3, 4, 5),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))
 is.validasterdata(fred)

 cmat <- constancy(fred, parm.type = "theta")
 multinomial.ind <- sort(unique(fred$redata$id[fred$recode == 1]))
 mycmat <- matrix(NA, nrow = length(multinomial.ind), ncol = nrow(fred$redata))
 for (ind in seq(along = multinomial.ind))
     mycmat[ind, ] <- as.numeric(fred$redata$id == ind & fred$recode == 1)
 identical(cmat, mycmat)
 
 sally <- fred
 y <- sally$redata[[sally$response.name]]
 v <- as.character(sally$redata$varb)
 ok <- y == 0 & v == "m1"
 i.save <- numeric(0)
 cmat.extra <- NULL
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     i.save <- i
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 y.buddy <- c(NA, y)[sally$regroup + 1]
 y.buddy.buddy <- c(NA, y.buddy)[sally$regroup + 1]
 ok <- y == 0 & y.buddy == 0 & v == "m2" & (! is.element(sally$regroup, i.save))
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     sally$redelta[sally$regroup[i]] <- (-1)
     i.save <- c(i.save, sally$regroup[i])
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
     foo <- rep(0, length(y))
     foo[sally$regroup[i]] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- seq(along = y)
 if (length(i.save) > 0)
    ok <- ok[- i.save]
 i <- ok[1]
 sally$redelta[c(i, i + 100, i + 200)] <- (- 1)
 ok <- y == 0 & v == "b1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- y == 1 & v == "b1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (+1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- y == 0 & v == "p1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- y == 1 & v == "z1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 
 cmat.sally <- constancy(sally, parm.type = "theta")
 cmat.char <- apply(cmat, 1, paste, collapse = "*")
 cmat.sally.char <- apply(cmat.sally, 1, paste, collapse = "*")
 cmat.extra.char <- apply(cmat.extra, 1, paste, collapse = "*")
 foo <- match(cmat, cmat.sally)
 all(! is.na(foo))
 identical(sort(cmat.sally), sort(c(cmat, cmat.extra)))

 ##### is.same parm.type = "theta"

 theta1 <- rnorm(length(fred$repred))
 theta2 <- rnorm(length(fred$repred))
 ! is.same(theta1, theta2, fred)
 theta2 <- theta1 + t(cmat) %*% rnorm(nrow(cmat))
 is.same(theta1, theta2, fred)

 ##### parm.type = "phi"

 rm(list = ls())

 data(echinacea)
 cmat <- constancy(echinacea, parm.type = "phi")
 nrow(cmat) == 0

 data(test1)
 fred <- asterdata(test1,
     vars = c("m1", "m2", "m3", "n1", "n2", "b1", "p1", "z1"),
     pred = c(0, 0, 0, 1, 1, 2, 3, 6), group = c(0, 1, 2, 0, 4, 0, 0, 0),
     code = c(1, 1, 1, 2, 2, 3, 4, 5),
     families = list(fam.multinomial(3), "normal.location.scale",
     "bernoulli", "poisson", "zero.truncated.poisson"))
 is.validasterdata(fred)

 cmat <- constancy(fred, parm.type = "phi")
 multinomial.ind <- sort(unique(fred$redata$id[fred$recode == 1]))
 mycmat <- matrix(NA, nrow = length(multinomial.ind), ncol = nrow(fred$redata))
 for (ind in seq(along = multinomial.ind))
     mycmat[ind, ] <- as.numeric(fred$redata$id == ind & fred$recode == 1)
 identical(cmat, mycmat)
 
 sally <- fred
 y <- sally$redata[[sally$response.name]]
 v <- as.character(sally$redata$varb)
 ok <- y == 0 & v == "m1"
 i.save <- numeric(0)
 cmat.extra <- NULL
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     i.save <- i
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 y.buddy <- c(NA, y)[sally$regroup + 1]
 y.buddy.buddy <- c(NA, y.buddy)[sally$regroup + 1]
 ok <- y == 0 & y.buddy == 0 & v == "m2" & (! is.element(sally$regroup, i.save))
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     sally$redelta[sally$regroup[i]] <- (-1)
     i.save <- c(i.save, sally$regroup[i])
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
     foo <- rep(0, length(y))
     foo[sally$regroup[i]] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- seq(along = y)
 if (length(i.save) > 0)
    ok <- ok[- i.save]
 i <- ok[1]
 sally$redelta[c(i, i + 100, i + 200)] <- (- 1)
 ok <- y == 0 & v == "b1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- y == 1 & v == "b1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (+1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     foo[sally$repred[i]] <- (-1)
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- y == 0 & v == "p1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     cmat.extra <- rbind(cmat.extra, foo)
 }
 ok <- y == 1 & v == "z1"
 if (any(ok)) {
     i <- seq(along = ok)[ok]
     i <- i[1]
     sally$redelta[i] <- (-1)
     foo <- rep(0, length(y))
     foo[i] <- 1
     foo[sally$repred[i]] <- (-1)
     cmat.extra <- rbind(cmat.extra, foo)
 }
 
 cmat.sally <- constancy(sally, parm.type = "phi")
 cmat.char <- apply(cmat, 1, paste, collapse = "*")
 cmat.sally.char <- apply(cmat.sally, 1, paste, collapse = "*")
 cmat.extra.char <- apply(cmat.extra, 1, paste, collapse = "*")
 foo <- match(cmat, cmat.sally)
 all(! is.na(foo))
 identical(sort(cmat.sally), sort(c(cmat, cmat.extra)))

 ##### is.same parm.type = "phi"

 phi1 <- rnorm(length(fred$repred))
 phi2 <- rnorm(length(fred$repred))
 ! is.same(phi1, phi2, fred, parm.type = "phi")
 phi2 <- phi1 + t(cmat) %*% rnorm(nrow(cmat))
 is.same(phi1, phi2, fred, parm.type = "phi")

