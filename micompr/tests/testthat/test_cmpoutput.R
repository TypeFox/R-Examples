# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

library(micompr)
context("cmpoutput")

test_that("cmpoutput constructs the expected objects", {

  # Minimum percentage of variance to be explained (< 1) OR
  # number of PCs specified directly (> 1)
  minvar <- 0.9
  ve_npcs <- c(0.3, 0.5, 0.7, 0.9, 3, 6, 9)

  # Which of these are variances?
  idxvar <- which(ve_npcs < 1)
  # And which are PCs?
  idxpcs <- which(ve_npcs > 1)

  # Instantiate several cmpoutput objects for testing using the datasets
  # provided with the package
  cmp_ok1 <- cmpoutput("SheepPop",
                       minvar,
                       pphpc_ok$data[["Pop.Sheep"]],
                       pphpc_ok$obs_lvls)

  cmp_noshuff2 <- cmpoutput("WolfPop",
                            minvar,
                            pphpc_noshuff$data[["Pop.Wolf"]],
                            pphpc_noshuff$obs_lvls,
                            lim_npcs = FALSE,
                            mnv_test = "Roy")

  cmp_diff7 <- cmpoutput("Concat",
                         minvar,
                         pphpc_diff$data[["All"]],
                         pphpc_diff$obs_lvls,
                         lim_npcs = FALSE)

  cmp_vlo6 <- cmpoutput("GrassEn",
                        minvar,
                        pphpc_testvlo$data[["Energy.Grass"]],
                        pphpc_testvlo$obs_lvls,
                        mnv_test = "Wilks")

  cmp_multve <- cmpoutput("WolfEn",
                          ve_npcs,
                          pphpc_ok$data[["Energy.Wolf"]],
                          pphpc_ok$obs_lvls,
                          mnv_test = "Hotelling-Lawley")

  # Instantiate a cmpoutput object with output from four pphpc implementations
  # Instantiate a grpobjects first
  go_quad <-
    grpoutputs(6,
               c(system.file("extdata", "nl_ok", package = "micompr"),
                 system.file("extdata", "j_ex_ok", package = "micompr"),
                 system.file("extdata", "j_ex_noshuff", package = "micompr"),
                 system.file("extdata", "j_ex_diff", package = "micompr")),
               rep(glob2rx("stats400v1*.tsv"), 4))

  cmp_quad3 <- cmpoutput("GrassQty",
                         minvar,
                         go_quad$data[["out3"]],
                         go_quad$obs_lvls)

  #### Start testing ####

  ## Common tests for the five cmpoutput objects ##
  for (ccmp in list(cmp_ok1, cmp_noshuff2, cmp_diff7,
                    cmp_vlo6, cmp_multve , cmp_quad3)) {

    # Check if cmpoutput objects have the correct type
    expect_is(ccmp, "cmpoutput")

    if (length(ccmp$ve) == 1) {

      # Test if the minimum percentage of variance to be explained matches what
      # was set at instantiation time
      expect_equal(ccmp$ve, minvar)

      # Check that the number of PCs which explain the specified minimum
      # percentage of variance has the expected value
      expect_equal(ccmp$npcs,
                   match(TRUE, cumsum(ccmp$varexp) > minvar))

    } else {

      # Test if the minimum percentage of variance to be explained OR if the
      # number of PCs matches what was set at instantiation time
      expect_equal(ccmp$ve[idxvar], ve_npcs[idxvar])
      expect_equal(ccmp$npcs[idxpcs], ve_npcs[idxpcs])

      # Check that the number of PCs which explain the specified minimum
      # percentage of variance has the expected value

      # Case 1 - Variance to explain was specified
      expect_equal(ccmp$npcs[idxvar],
                   sapply(ve_npcs[idxvar],
                          function(mv, ve) match(T, cumsum(ve) > mv),
                          ccmp$varexp))
      # Case 2 - Number of PCs was directly specified
      expect_equal(ccmp$ve[idxpcs],
                   sapply(ve_npcs[idxpcs],
                          function(npcs, ve) sum(ve[1:npcs]),
                          ccmp$varexp))

    }

    # Check that the tests objects are what is expected and p-values are within
    # the 0-1 range
    for (i in 1:length(ccmp$npcs)) {

      if (ccmp$npcs[i] > 1) {
        expect_is(ccmp$tests$manova[[i]], "manova")
        expect_true((ccmp$p.values$manova[i] >= 0)
                    && (ccmp$p.values$manova[i] <= 1),
                    "MANOVA p-value not between 0 and 1.")
      } else {
        if (length(ccmp$ccmp$tests$manova) > 1) {
          expect_null(ccmp$tests$manova[[i]])
        }
      }
    }

    for (i in length(ccmp$varexp)) {
      if (length(levels(ccmp$obs_lvls)) == 2) {
        expect_is(ccmp$tests$parametric[[i]], "htest")
      } else {
        expect_is(ccmp$tests$parametric[[i]], "aov")
      }
      expect_is(ccmp$tests$nonparametric[[i]], "htest")
    }
    expect_true((ccmp$p.values$parametric[i] >= 0)
                && (ccmp$p.values$parametric[i] <= 1),
                "Parametric test p-value not between 0 and 1.")
    expect_true((ccmp$p.values$nonparametric[i] >= 0)
                && (ccmp$p.values$nonparametric[i] <= 1),
                "Non-parametric test p-value not between 0 and 1.")
    expect_equal(ccmp$p.values$parametric_adjusted,
                 pmin(ccmp$p.values$parametric / ccmp$varexp, 1))
    expect_equal(ccmp$p.values$nonparametric_adjusted,
                 pmin(ccmp$p.values$nonparametric / ccmp$varexp, 1))

  }
  ## Different tests for the cmpoutput objects ##

  # Check the names given to the comparisons
  expect_equal(cmp_ok1$name, "SheepPop")
  expect_equal(cmp_noshuff2$name, "WolfPop")
  expect_equal(cmp_diff7$name, "Concat")
  expect_equal(cmp_vlo6$name, "GrassEn")
  expect_equal(cmp_quad3$name, "GrassQty")

  # Check that the dimensions of the scores (PCA) matrix are limited by the
  # number of observations
  # In these cases output length (number of variables) is always larger than
  # the number of observations
  expect_equal(dim(cmp_ok1$scores),
               c(length(pphpc_ok$obs_lvls), length(pphpc_ok$obs_lvls)))
  expect_equal(dim(cmp_noshuff2$scores),
               c(length(pphpc_noshuff$obs_lvls),
                 length(pphpc_noshuff$obs_lvls)))
  expect_equal(dim(cmp_diff7$scores),
               c(length(pphpc_diff$obs_lvls), length(pphpc_diff$obs_lvls)))
  expect_equal(dim(cmp_vlo6$scores),
               c(length(pphpc_testvlo$obs_lvls), length(pphpc_testvlo$obs_lvls)))
  expect_equal(dim(cmp_quad3$scores),
               c(length(go_quad$obs_lvls), length(go_quad$obs_lvls)))

  # Check if the observation levels are the same as in the original data
  expect_equal(cmp_ok1$obs_lvls, pphpc_ok$obs_lvls)
  expect_equal(cmp_noshuff2$obs_lvls, pphpc_noshuff$obs_lvls)
  expect_equal(cmp_diff7$obs_lvls, pphpc_diff$obs_lvls)
  expect_equal(cmp_vlo6$obs_lvls, pphpc_testvlo$obs_lvls)
  expect_equal(cmp_quad3$obs_lvls, go_quad$obs_lvls)

})


test_that("cmpoutput throws errors when improperly invoked", {

  # Test for incorrect 've' parameter
  expect_error(
    cmpoutput("B", -0.01, pphpc_ok$data[[2]], pphpc_ok$obs_lvls),
    "'ve_npcs' parameter must only have positive values.",
    fixed = TRUE
  )

  # Test for invalid number of levels
  expect_error(
    cmpoutput("C", 0.5, pphpc_ok$data[[3]], rep(1, length(pphpc_ok$obs_lvls))),
    "At least two levels are required to perform model comparison.",
    fixed = TRUE
  )

  # Test for error due to different number of observations in data and
  # observation levels
  expect_error(
    cmpoutput("D", 0.3, pphpc_ok$data[[4]], rep(pphpc_ok$obs_lvls, 2)),
    "Number of observations in 'data' and 'obs_lvls' does not match.",
    fixed = TRUE
  )

})

test_that("assumptions.cmpoutput creates the correct object", {

  #### No warnings #####

  # Create a cmpoutput object from the provided datasets
  cmp <- cmpoutput("All", 0.8, pphpc_ok$data[["All"]], pphpc_ok$obs_lvls)

  # Get the assumptions for the parametric tests performed in cmp
  acmp <- assumptions(cmp)

  # Check that the objects are of the correct type
  expect_is(acmp, "assumptions_cmpoutput")
  for (amnv in acmp$manova) {
    if (!is.null(amnv)) {
      expect_is(amnv, "assumptions_manova")
    }
  }
  expect_is(acmp$ttest, "assumptions_paruv")

  #### Warnings about more variables than observations ####

  # Create a cmpoutput object from the provided datasets, set ve to 0.9 such
  # that more PCs are required than before
  cmp <- cmpoutput("All", 0.9, pphpc_ok$data[["All"]], pphpc_ok$obs_lvls)

  # Check that the creation of the assumptions object produces the expected
  # warnings, i.e. that there are more variables (represented by the PCs) than
  # observations
  expect_warning(assumptions(cmp),
                 paste("Royston test requires more observations than",
                       "(dependent) variables (DVs). Reducing number of",
                       "variables from 10 to 9 in group"),
                 fixed = TRUE)

  # Get the assumptions for the parametric tests performed in cmp, disable
  # warnings this time
  oldw <- getOption("warn")
  options(warn = -1)

  acmp <- assumptions(cmp)

  options(warn = oldw)

  # Check that the objects are of the correct type
  expect_is(acmp, "assumptions_cmpoutput")
  for (amnv in acmp$manova) {
    if (!is.null(amnv)) {
      expect_is(amnv, "assumptions_manova")
    }
  }
  expect_is(acmp$ttest, "assumptions_paruv")

  #### Test with insufficient observations for the Royston test ####

  # In this case it should not be possible to perform the Royston test
  cmp <- cmpoutput("GrassEn",
                   0.9,
                   pphpc_testvlo$data[["Energy.Grass"]],
                   pphpc_testvlo$obs_lvls)

  # Check that if it is so
  expect_warning(assumptions(cmp),
                 "Royston test requires at least 4 observations",
                 fixed = TRUE)

  # Do create the object and check remaining stuff
  oldw <- getOption("warn")
  options(warn = -1)

  acmp <- assumptions(cmp)

  options(warn = oldw)

  # Check that the objects are of the correct type
  expect_is(acmp, "assumptions_cmpoutput")
  for (amnv in acmp$manova) {
    if (!is.null(amnv)) {
      expect_is(amnv, "assumptions_manova")
    }
  }
  expect_is(acmp$ttest, "assumptions_paruv")


})