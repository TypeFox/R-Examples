library("geiger")
geo <- get(data(geospiza))
tmp <- treedata(geo$phy, geo$dat)
phy <- tmp$phy
dat <- tmp$data[,1]
nboot <- 10
mc.cores <- 1




test_that("we can run pmc_fit", {
  fit_A <- pmc_fit(tree = phy, data = dat, model = "BM")
  expect_is(fit_A, "gfit")
})

test_that("we can run tidy_gather internal methods", {
  fit_A <- pmc_fit(tree = phy, data = dat, model = "BM")
  A_sims <- format_sims(simulate(fit_A, nboot))
  
  AA <- parallel::mclapply(1:nboot, function(i) update(fit_A, A_sims[,i]), mc.cores = mc.cores)
  
  mtrx <- sapply(AA, function(x) {
    out <- coef(x)
    if(is.list(out))
      out <- unlist(out)
    out
  })
  tmp <- data.frame(t(rbind(mtrx, rep = 1:dim(mtrx)[[2]])))
  who <- which(names(tmp)!="rep")
  
  tidy_AA <- data.frame(comparison = "AA", tidyr::gather_(tmp, "parameter", "value", names(tmp)[who]))
  expect_is(tidy_AA, "data.frame")
})



test_that("we can run tidy_gather", {
  fit_A <- pmc_fit(tree = phy, data = dat, model = "BM")
  A_sims <- format_sims(simulate(fit_A, nboot))
  
  AA <- parallel::mclapply(1:nboot, function(i) update(fit_A, A_sims[,i]), mc.cores = mc.cores)
  tidy_AA <- tidy_pars(AA)
  expect_is(tidy_AA, "data.frame")
})

test_that("run pmc",
          {
            
            out <- pmc(phy, dat, "BM", "lambda", nboot = 10, mc.cores = 1)
            expect_is(out, "pmc")

            p <- plot(out)
            expect_is(p, "ggplot") 
          })
