context("wjd")

test_that("wjd() gives results similar to White et al., 1958", {
  # the values from last column of Table III in the paper
  X <- c(0.040668, 0.147730, 0.783153, 0.001414, 0.485247, 0.000693, 0.027399, 0.017947, 0.037314, 0.096872)
  w <- wjd()
  expect_equal(X, w$X, tolerance=1e-4)
})

test_that("guess() operates on intermediate compositions but fails with endmembers", {
  alkanes <- c("n-hexane", "n-heptane", "n-octane", "n-nonane")
  ialk <- info(alkanes, "liq")
  expect_error(run.guess(ialk, "C6H14"), "there are only 0")
  # hmm, on windows this has a length of 5 (20120626)
  # probably should filter out guesses with very low abundances
  #expect_true(length(run.guess(ialk, "C7H16"))==4)
  expect_true(length(run.guess(ialk, "C8H18"))==5)
  expect_error(run.guess(ialk, "C9H20"), "there are only 0")
})

test_that("open-system equilibrium distributions reproduce the results of wjd()", {
  ### set up system
  # use proteins in the lipid particle (n=19)
  y <- yeastgfp("lipid.particle")
  # get the amino acid compositions of the proteins
  aa <- more.aa(y$protein, "Sce")
  # don't use those with NA abundance or sequence (leaves n=17)
  ina <- is.na(y$abundance) | is.na(aa$chains)
  aa <- aa[!ina, ]
  # normalize the proteins to single residues; columns 6:25 are the amino acid counts
  aa.625 <- aa[, 6:25]
  aa[, 6:25] <- aa.625 / rowSums(aa.625)
  # add proteins to thermo$protein
  iprotein <- add.protein(aa)
  # add proteins to thermo$obigt
  iobigt <- info(paste(aa$protein, aa$organism, sep="_"))
  ### closed system calculation (constant composition of elements)
  # use equal initial abundances
  Y <- rep(100, length(iobigt))
  # run the Gibbs energy minimization (this did not iterate before 20130109,
  # due to bug in calculation of free energy derivative)
  w <- run.wjd(iobigt, Y=Y, Gfrac=1e-15, nlambda=1001)
  # the molar abundances
  X.closed <- w$X
  # get the chemical potentials of the elements
  ep <- equil.potentials(w)
  # the corresponding logarithms of activities of the basis species
  basis("CHNOS")
  bl <- basis.logact(ep)
  ### open system calculation (constant chemical potentials of basis species)
  # set the activities of the basis species
  basis(names(bl), bl)
  # get the affinities of the formation reactions
  a <- affinity(iprotein=iprotein)
  # then the equilibrium abundances, with total moles of residues as used in wjd
  e <- equilibrate(a, loga.balance=log10(sum(Y)))
  X.open <- 10^unlist(e$loga.equil)
  # the test: abundances calculated both ways are equal
  expect_equal(X.closed, X.open, tolerance=0.019)
  # seems that we could do better than that 1.9% mean difference!
})

# see also test-swap.basis.R for an example using run.wjd() and 
# equil.potentials() to generate chemical potentials of elements

# references

# White, W. B., Johnson, S. M. and Dantzig, G. B. (1958) 
#   Chemical equilibrium in complex mixtures. 
#   J. Chem. Phys. 28, 751--755. http://dx.doi.org/10.1063/1.1744264
