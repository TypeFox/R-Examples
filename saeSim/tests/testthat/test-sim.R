context("sim")
test_that("Method for base", {
  dat <- sim_base(base_id(nDomains=3, nUnits = 4)) %>%
    # Fixed-Effects: drawn from N(50, 20^2), b0 = 0, slope = 10
    sim_gen_x(mean=50, sd=20, name = "x") %>%
    # Model-Error: e ~ N(0, 1)
    sim_gen_e(0, 1) %>%
    # Random-Intercept: v ~ N(0, 1)
    sim_gen_v(0, 1) %>%
    # Correlated random-effect following a SAR(1)
    sim_gen(gen_v_sar(mean=0, sd=2, rho=0.5, type="rook", name = "v_sp")) %>%
    # Adding outliers in each error-component:
    # 5% in each area,at least 1, those with highest ID:
    sim_gen_ec(sd = 150, nCont=0.05, type="unit", areaVar = "idD", fixed=TRUE) %>%
    # 5% of the areas, at least 1, those with highest ID
    sim_gen_vc(sd = 50, nCont=0.05, type="area", areaVar = "idD", fixed=TRUE) %>%
    # 2 Areas, randomly chosen
    sim_gen_cont(gen_v_sar(sd = 50, name = "v_sp"), nCont=2, type="area", areaVar = "idD", fixed=FALSE) %>%
    # 1 in area1, 2 in area2, 1 in area3
    sim_gen_vc(mean = 10, sd = 1, nCont = c(1, 2, 1), type = "unit", areaVar = "idD", fixed = TRUE) %>%
    # 2 outliers, somewhere...
    sim_gen_vc(mean = 10, sd = 1, nCont = 2, type = "unit", areaVar = NULL, fixed = FALSE) %>%
    as.data.frame
  
  expect_equal(nrow(dat), 12)
  expect_equal(ncol(dat), 7)
  expect_true(all(c("idU", "idD") %in% names(dat)))
  
})

test_that("Method for sim_setup", {
  setup <- sim_base(base_id(nDomains = 4, nUnits = 3)) %>%
    sim_gen_x() %>% sim_gen_e()
  datList <- sim(setup %>% sim_simName("test"), R = 500)
  
  expect_equal(length(datList), (500))
  expect_equal(max(rbind_all(datList)$idR), (500))
  expect_that(all(rbind_all(datList)$simName == "test"), is_true())
  
})

test_that("Save incremental files", {
  
  tmp <- tempdir()
  sim_base_lm() %>%
    sim(R = 2, path = tmp)
  
  expect_equal(nrow(sim_read_data(tmp)), 20000)
  expect_equal(ncol(sim_read_data(tmp)), 7)
  
  sim_base_lm() %>%
    sim_agg(function(dat) list(dat = dat, el = 1)) %>%
    sim(R = 2, path = tmp, fileExt = ".RData")
  
  expect_equal(nrow(sim_read_list(tmp)[[1]]$dat), 10000)
  expect_equal(ncol(sim_read_list(tmp)[[1]]$dat), 5)
  expect_equal(length(sim_read_list(tmp)[[1]]), 4)
  expect_equal(length(sim_read_list(tmp)), 2)
})
