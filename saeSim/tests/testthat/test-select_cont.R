context("select_cont")
test_that(desc = "select_cont", {
  dat1 <- sim_base(base_id(nDomains=5, nUnits = 2)) %>% sim_gen_e() %>% as.data.frame
  expect_equal(sum(select_cont(dat1, 1L, "unit", "idD", TRUE)$e == 0), 5)
  expect_equal(sum(select_cont(dat1, 1L, "unit", "idD", FALSE)$e == 0), 5)
  expect_equal(sum(select_cont(dat1, 1L, "area", "idD", TRUE)$e == 0), 8)
  expect_equal(sum(select_cont(dat1, 1L, "area", "idD", FALSE)$e == 0), 8)
  expect_equal(sum(select_cont(dat1, 1L, "unit", NULL, TRUE)$e == 0), 9)
  expect_equal(sum(select_cont(dat1, 1L, "unit", NULL, FALSE)$e == 0), 9)
  
  expect_true(sum(select_cont(dat1, 0.5, "unit", "idD", TRUE)$e == 0) > 0)
  expect_true(sum(select_cont(dat1, 0.5, "area", "idD", TRUE)$e == 0) > 0)
  expect_true(sum(select_cont(dat1, 0.5, "unit", NULL, TRUE)$e == 0) > 0)
  
  dat2 <- sim_base(base_id(nDomains = 5, nUnits = 1:5)) %>% sim_gen_e() %>% as.data.frame
  expect_equal(sum(select_cont(dat2, 1L, "unit", "idD", TRUE)$e == 0), 10)
  expect_equal(sum(select_cont(dat2, 1L, "unit", "idD", FALSE)$e == 0), 10)
  expect_equal(sum(select_cont(dat2, 1L, "area", "idD", TRUE)$e == 0), 10)
  expect_equal(sum(select_cont(dat2, 1L, "unit", NULL, TRUE)$e == 0), 14)
  expect_equal(sum(select_cont(dat2, 1L, "unit", NULL, FALSE)$e == 0), 14)
  
  numberOfContObs <- (select_cont(dat2, nCont = c(0, 0, 0, 0.2, 1), "unit", "idD", TRUE)$idC %>% sum)
  (numberOfContObs >= 5 & numberOfContObs <= 9) %>%
    expect_true
  
  testthat::expect_equal(sum(select_cont(dat2, c(1, 4), "area", "idD", TRUE)$e == 0), 10)
  
})
