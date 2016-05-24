context("Replicating results from Behavior Research Manuscript")

test_that("Experiment 1", {
  # load data
  data(exp1)  
  # this uses code from likelihood_d() to calculate the mean complexity K
  # for all strings of length 10 with alphabet = 4:
  tmp <- acss_data[nchar(rownames(acss_data)) == 10, "K.4", drop = FALSE]
  tmp <- tmp[!is.na(tmp[,"K.4"]),,drop = FALSE]
  tmp$count <- count_class(rownames(tmp), alphabet = 4)
  mean_K <- with(tmp, sum(K.4*count)/sum(count))  
  res_t <- t.test(acss(exp1$string, 4)[,"K.4"], mu = mean_K)
  expect_that(res_t$p.value, is_less_than(0.001))
  expect_that(round(res_t$statistic, 2), is_equivalent_to(10.62))
})

test_that("Experiment 2", {
  data(exp2)
  exp2$K <- acss(exp2$string, 6)[,"K.6"]
  m_log <- glm(response ~ K, exp2, family = binomial)
  # slope
  expect_equivalent(round(coef(m_log)[2], 1), 1.9)
  # odds ratio of K:
  expect_equivalent(round(exp(coef(m_log)[2]), 2), 6.69)  
  # subjective threshold
  expect_equivalent(round(-coef(m_log)[1]/coef(m_log)[2], 2), 36.74)
  rm(exp2)
})

test_that("Matthews 2013", {
  data(matthews2013)
  spans <- c(4, 5, 10)
  # note, the next loop takes more than 5 minutes.
  for (i in spans) {
    matthews2013[,paste0("K2_span", i)] <- 
      sapply(local_complexity(matthews2013$string, alphabet=2, span = i), mean)
  }
  
  lm_list <- vector("list", 8)
  for (i in seq_along(spans)) {
    lm_list[[i]] <- lm(as.formula(paste0("mean ~ K2_span", spans[i])), matthews2013)
  }
  
  expect_equivalent(round(summary(lm_list[[1]])$r.squared, 2), 0.54)
  expect_equivalent(round(summary(lm_list[[2]])$r.squared, 2), 0.50)
  expect_less_than(summary(lm_list[[3]])$r.squared, 0.00002)
  rm(matthews2013)
})