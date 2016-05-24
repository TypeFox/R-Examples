
 library(aster)

 data(radish)

 options(digits=4) # avoid rounding differences

 pred <- c(0,1,2)
 fam <- c(1,3,2)

 ### need object of type aster to supply to penmlogl and pickle

 aout <- aster(resp ~ varb + fit : (Site * Region),
     pred, fam, varb, id, root, data = radish)

 ### model matrices for fixed and random effects

 modmat.fix <- model.matrix(resp ~ varb + fit : (Site * Region), data = radish)
 modmat.blk <- model.matrix(resp ~ 0 + fit:Block, data = radish)
 modmat.pop <- model.matrix(resp ~ 0 + fit:Pop, data = radish)

 rownames(modmat.fix) <- NULL
 rownames(modmat.blk) <- NULL
 rownames(modmat.pop) <- NULL

 idrop <- match(aout$dropped, colnames(modmat.fix))
 idrop <- idrop[! is.na(idrop)]
 modmat.fix <- modmat.fix[ , - idrop]

 nfix <- ncol(modmat.fix)
 nblk <- ncol(modmat.blk)
 npop <- ncol(modmat.pop)

 ### try reaster

 rout <- reaster(resp ~ varb + fit : (Site * Region),
     list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
     pred, fam, varb, id, root, data = radish)
 summary(rout)
 summary(rout, stand = FALSE)

 class(rout)
 names(rout)

 identical(modmat.fix, rout$fixed)
 identical(list(block = modmat.blk, pop = modmat.pop), rout$random)
 class(aout)
 class(rout$obj)
 foo <- rep(NA, length(rout$obj))
 names(foo) <- names(rout$obj)
 for (n in names(rout$obj))
     foo[n] <- identical(aout[[n]], rout$obj[[n]])
 foo

 foo <- aout$x
 is.matrix(foo)
 dimnames(foo) <- NULL
 identical(foo, rout$obj$x)
 foo <- aout$root
 is.matrix(foo)
 dimnames(foo) <- NULL
 identical(foo, rout$obj$root)
 foo <- aout$modmat
 is.array(foo)
 dimnames(foo) <- list(NULL, NULL, dimnames(foo)[[3]])
 identical(foo, rout$obj$modmat)

 alpha.mle <- rout$alpha
 sigma.mle <- rout$sigma
 c.mle <- rout$c
 fixed <- rout$fixed
 random <- rout$random
 bee.mle <- rout$b
 nu.mle <- rout$nu

 zwz.mle <- makezwz(sigma.mle, c(alpha.mle, c.mle), fixed = fixed,
     random = random, obj = rout$obj)
 identical(zwz.mle, rout$zwz)

 qout <- quickle(c(alpha.mle, nu.mle), bee.mle, fixed = fixed,
     random = random, obj = rout$obj, zwz = zwz.mle, deriv = 2)
 sout <- summary(rout)
 identical(qout$hessian, sout$fish)

 # now check standard error calculations in summary.reaster

 fred <- eigen(sout$fish, symmetric = TRUE)
 fish.inv <- fred$vectors %*% diag(1 / fred$values) %*% t(fred$vectors)
 all.equal(fish.inv, solve(sout$fish))
 se.parm <- sqrt(diag(fish.inv))
 se.alpha <- se.parm[seq(along = se.parm) <= nfix]
 se.nu <- se.parm[seq(along = se.parm) > nfix]
 identical(se.alpha, as.vector(sout$alpha[ , 2]))
 identical(se.nu, as.vector(sout$nu[ , 2]))
 identical(se.nu / (2 * sigma.mle), sout$sigma[ , 2])

 eps <- 1e-4
 npoin <- 9
 stopifnot(npoin %% 2 == 1)
 alphanu.mle <- c(alpha.mle, nu.mle)
 bsp <- matrix(NA_real_, nrow = length(bee.mle), ncol = length(alphanu.mle))
 for (i in seq(along = alphanu.mle)) {
     blurfle <- matrix(NA_real_, nrow = length(bee.mle), ncol = npoin)
     foo <- alphanu.mle
     bar <- double(npoin)
     for (j in 1:npoin) {
         foo[i] <- alphanu.mle[i] + (j - (npoin + 1) / 2) * eps
         bar[j] <- foo[i]
         qout <- quickle(foo, bee.mle, fixed = fixed,
             random = random, obj = rout$obj, zwz = zwz.mle)
         blurfle[ , j] <- qout$bee
     }
     for (j in seq(along = bee.mle)) {
         surfle <- splinefun(bar, blurfle[j, ], method = "natural")
         bsp[j, i] <- surfle(alphanu.mle[i], deriv = 1)
     }
 }

 qout <- quickle(c(alpha.mle, nu.mle), bee.mle, fixed = fixed,
     random = random, obj = rout$obj, zwz = zwz.mle, deriv = 2)
 bsp2 <- with(qout, pbb.inv %*% cbind(pba, pbn))
 all.equal(- bsp, bsp2, tol = 5e-6)

 ### try reaster with no data frame

 y <- radish$resp
 v <- radish$varb
 f <- radish$fit
 s <- radish$Site
 r <- radish$Region
 b <- radish$Block
 p <- radish$Pop
 i <- radish$id
 rt <- radish$root
 rout <- reaster(y ~ v + f : (s * r),
     list(block = ~ 0 + f : b, pop = ~ 0 + f : p),
     pred, fam, v, i, rt)
 summary(rout)

