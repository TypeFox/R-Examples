context("tableBalance")

################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 13.11.2015: Updated hb=2 to hb=3 in to correspond to new implemented method in calculateBalance.
# 07.05.2014: First version.
# 
# 
# test_dir("inst/tests/")
# test_file("tests/testthat/test-tableBalance.r")
# test_dir("tests/testthat")

test_that("tableBalance", {

  # Get test data.
  data(set2)
  data(ref2)

  # TEST 01 -------------------------------------------------------------------
  
  # Analyse dataframe.
  tmp <- calculateBalance(data=set2, ref=ref2, lb="prop",
                          per.dye=TRUE, hb=3,
                          ignore.case=TRUE)
  
  res <- tableBalance(data=tmp, scope="locus", quant=0.05)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$Hb.n))
  expect_false(is.null(res$Hb.Mean))
  expect_false(is.null(res$Hb.Stdv))
  expect_false(is.null(res$Hb.Perc.5))
  expect_false(is.null(res$Lb.n))
  expect_false(is.null(res$Lb.Mean))
  expect_false(is.null(res$Lb.Stdv))
  expect_false(is.null(res$Lb.Perc.5))
  
  # Check for NA's.
  expect_false(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Hb.n)))
  expect_true(any(is.na(res$Hb.Min)))
  expect_true(any(is.na(res$Hb.Mean)))
  expect_true(any(is.na(res$Hb.Stdv)))
  expect_true(any(is.na(res$Hb.Perc.5)))
  expect_false(any(is.na(res$Lb.n)))
  expect_false(any(is.na(res$Lb.Min)))
  expect_false(any(is.na(res$Lb.Mean)))
  expect_true(any(is.na(res$Lb.Stdv)))
  expect_false(any(is.na(res$Lb.Perc.5)))
  
  # Check result: number of values for heterozygous balance.
  expect_that(res$Hb.n[1], equals(2))
  expect_that(res$Hb.n[2], equals(0))
  expect_that(res$Hb.n[3], equals(2))
  expect_that(res$Hb.n[4], equals(0))
  expect_that(res$Hb.n[5], equals(2))
  expect_that(res$Hb.n[6], equals(0))
  expect_that(res$Hb.n[7], equals(0))
  expect_that(res$Hb.n[8], equals(2))
  expect_that(res$Hb.n[9], equals(0))
  expect_that(res$Hb.n[10], equals(2))
  expect_that(res$Hb.n[11], equals(0))

  # Check result: Mean heterozygous balance.
  expect_that(res$Hb.Mean[1], equals(mean(c(402/460,215/225))))
  expect_that(res$Hb.Mean[2], equals(as.numeric(NA)))
  expect_that(res$Hb.Mean[3], equals(mean(c(423/491,198/241))))
  expect_that(res$Hb.Mean[4], equals(as.numeric(NA)))
  expect_that(res$Hb.Mean[5], equals(mean(c(587/632,312/326))))
  expect_that(res$Hb.Mean[6], equals(as.numeric(NA)))
  expect_that(res$Hb.Mean[7], equals(as.numeric(NA)))
  expect_that(res$Hb.Mean[8], equals(mean(c(361/398,195/206))))
  expect_that(res$Hb.Mean[9], equals(as.numeric(NA)))
  expect_that(res$Hb.Mean[10], equals(mean(c(359/384,179/183))))
  expect_that(res$Hb.Mean[11], equals(as.numeric(NA)))

  # Check result: Heterozygous balance standard deviation.
  expect_that(res$Hb.Stdv[1], equals(sd(c(402/460,215/225))))
  expect_that(res$Hb.Stdv[2], equals(as.numeric(NA)))
  expect_that(res$Hb.Stdv[3], equals(sd(c(423/491,198/241))))
  expect_that(res$Hb.Stdv[4], equals(as.numeric(NA)))
  expect_that(res$Hb.Stdv[5], equals(sd(c(587/632,312/326))))
  expect_that(res$Hb.Stdv[6], equals(as.numeric(NA)))
  expect_that(res$Hb.Stdv[7], equals(as.numeric(NA)))
  expect_that(res$Hb.Stdv[8], equals(sd(c(361/398,195/206))))
  expect_that(res$Hb.Stdv[9], equals(as.numeric(NA)))
  expect_that(res$Hb.Stdv[10], equals(sd(c(359/384,179/183))))
  expect_that(res$Hb.Stdv[11], equals(as.numeric(NA)))
  
  # Check result: 5 percentile for heterozygous balance.
  expect_that(res$Hb.Perc.5[1], equals(as.numeric(quantile(c(402/460,215/225), 0.05))))
  expect_that(res$Hb.Perc.5[2], equals(as.numeric(NA)))
  expect_that(res$Hb.Perc.5[3], equals(as.numeric(quantile(c(423/491,198/241), 0.05))))
  expect_that(res$Hb.Perc.5[4], equals(as.numeric(NA)))
  expect_that(res$Hb.Perc.5[5], equals(as.numeric(quantile(c(587/632,312/326), 0.05))))
  expect_that(res$Hb.Perc.5[6], equals(as.numeric(NA)))
  expect_that(res$Hb.Perc.5[7], equals(as.numeric(NA)))
  expect_that(res$Hb.Perc.5[8], equals(as.numeric(quantile(c(361/398,195/206), 0.05))))
  expect_that(res$Hb.Perc.5[9], equals(as.numeric(NA)))
  expect_that(res$Hb.Perc.5[10], equals(as.numeric(quantile(c(359/384,179/183), 0.05))))
  expect_that(res$Hb.Perc.5[11], equals(as.numeric(NA)))
  
  # Check result: number of values for locus balance.
  expect_that(res$Lb.n[1], equals(2))
  expect_that(res$Lb.n[2], equals(2))
  expect_that(res$Lb.n[3], equals(2))
  expect_that(res$Lb.n[4], equals(2))
  expect_that(res$Lb.n[5], equals(2))
  expect_that(res$Lb.n[6], equals(2))
  expect_that(res$Lb.n[7], equals(2))
  expect_that(res$Lb.n[8], equals(2))
  expect_that(res$Lb.n[9], equals(1))
  expect_that(res$Lb.n[10], equals(1))
  expect_that(res$Lb.n[11], equals(1))
  
  # Check result: Mean locus balance.
  expect_that(res$Lb.Mean[1], equals(mean(c(862/3005,440/1530))))
  expect_that(res$Lb.Mean[2], equals(mean(c(506/3005,304/1530))))
  expect_that(res$Lb.Mean[3], equals(mean(c(914/3005,439/1530))))
  expect_that(res$Lb.Mean[4], equals(mean(c(723/3005,347/1530))))
  expect_that(res$Lb.Mean[5], equals(mean(c(1219/3363,638/1750))))
  expect_that(res$Lb.Mean[6], equals(mean(c(619/3363,309/1750))))
  expect_that(res$Lb.Mean[7], equals(mean(c(766/3363,402/1750))))
  expect_that(res$Lb.Mean[8], equals(mean(c(759/3363,401/1750))))
  expect_that(res$Lb.Mean[9], equals(592/1760))
  expect_that(res$Lb.Mean[10], equals(743/1760))
  expect_that(res$Lb.Mean[11], equals(425/1760))

  
  # Check result: Heterozygous balance standard deviation.
  expect_that(res$Lb.Stdv[1], equals(sd(c(862/3005,440/1530))))
  expect_that(res$Lb.Stdv[2], equals(sd(c(506/3005,304/1530))))
  expect_that(res$Lb.Stdv[3], equals(sd(c(914/3005,439/1530))))
  expect_that(res$Lb.Stdv[4], equals(sd(c(723/3005,347/1530))))
  expect_that(res$Lb.Stdv[5], equals(sd(c(1219/3363,638/1750))))
  expect_that(res$Lb.Stdv[6], equals(sd(c(619/3363,309/1750))))
  expect_that(res$Lb.Stdv[7], equals(sd(c(766/3363,402/1750))))
  expect_that(res$Lb.Stdv[8], equals(sd(c(759/3363,401/1750))))
  expect_that(res$Lb.Stdv[9], equals(as.numeric(NA)))
  expect_that(res$Lb.Stdv[10], equals(as.numeric(NA)))
  expect_that(res$Lb.Stdv[11], equals(as.numeric(NA)))
  
  # Check result: 5 percentile for heterozygous balance.
  expect_that(res$Lb.Perc.5[1], equals(as.numeric(quantile(c(862/3005,440/1530), 0.05))))
  expect_that(res$Lb.Perc.5[2], equals(as.numeric(quantile(c(506/3005,304/1530), 0.05))))
  expect_that(res$Lb.Perc.5[3], equals(as.numeric(quantile(c(914/3005,439/1530), 0.05))))
  expect_that(res$Lb.Perc.5[4], equals(as.numeric(quantile(c(723/3005,347/1530), 0.05))))
  expect_that(res$Lb.Perc.5[5], equals(as.numeric(quantile(c(1219/3363,638/1750), 0.05))))
  expect_that(res$Lb.Perc.5[6], equals(as.numeric(quantile(c(619/3363,309/1750), 0.05))))
  expect_that(res$Lb.Perc.5[7], equals(as.numeric(quantile(c(766/3363,402/1750), 0.05))))
  expect_that(res$Lb.Perc.5[8], equals(as.numeric(quantile(c(759/3363,401/1750), 0.05))))
  expect_that(res$Lb.Perc.5[9], equals(as.numeric(quantile(c(592/1760), 0.05))))
  expect_that(res$Lb.Perc.5[10], equals(as.numeric(quantile(c(743/1760), 0.05))))
  expect_that(res$Lb.Perc.5[11], equals(as.numeric(quantile(c(425/1760), 0.05))))
  
  
  # TEST 02 -------------------------------------------------------------------
  
  # Analyse dataframe.
  tmp <- calculateBalance(data=set2, ref=ref2, lb="prop",
                          per.dye=TRUE, hb=3,
                          ignore.case=TRUE)
  
  res <- tableBalance(data=tmp, scope="global", quant=0.10)
  
  # Check return class.  
  expect_that(class(res), matches(class(data.frame())))
  
  # Check that expected columns exist.  
  expect_true(is.null(res$Sample.Name))
  expect_false(is.null(res$Marker))
  expect_false(is.null(res$Hb.n))
  expect_false(is.null(res$Hb.Mean))
  expect_false(is.null(res$Hb.Stdv))
  expect_false(is.null(res$Hb.Perc.10))
  expect_false(is.null(res$Lb.n))
  expect_false(is.null(res$Lb.Mean))
  expect_false(is.null(res$Lb.Stdv))
  expect_false(is.null(res$Lb.Perc.10))
  
  # Check for NA's.
  expect_true(any(is.na(res$Marker)))
  expect_false(any(is.na(res$Hb.n)))
  expect_false(any(is.na(res$Hb.Min)))
  expect_false(any(is.na(res$Hb.Mean)))
  expect_false(any(is.na(res$Hb.Stdv)))
  expect_false(any(is.na(res$Hb.Perc.10)))
  expect_false(any(is.na(res$Lb.n)))
  expect_false(any(is.na(res$Lb.Min)))
  expect_false(any(is.na(res$Lb.Mean)))
  expect_false(any(is.na(res$Lb.Stdv)))
  expect_false(any(is.na(res$Lb.Perc.10)))
  
  # Check result: number of values for heterozygous balance.
  expect_that(res$Hb.n[1], equals(10))
  
  # Check result: Mean heterozygous balance.
  expect_that(res$Hb.Mean[1], equals(mean(c(402/460,215/225,
                                            423/491,198/241, 
                                            587/632,312/326,
                                            361/398,195/206,
                                            359/384,179/183))))
  
  # Check result: Heterozygous balance standard deviation.
  expect_that(res$Hb.Stdv[1], equals(sd(c(402/460,215/225,
                                          423/491,198/241,
                                          587/632,312/326,
                                          361/398,195/206,
                                          359/384,179/183))))
  
  # Check result: 10 percentile for heterozygous balance.
  expect_that(res$Hb.Perc.10[1], equals(as.numeric(quantile(c(402/460,215/225,
                                                              423/491,198/241,
                                                              587/632,312/326,
                                                              361/398,195/206,
                                                              359/384,179/183), 0.10))))

  
  # Check result: number of values for locus balance.
  expect_that(res$Lb.n[1], equals(19))
  
  # Check result: Mean locus balance.
  expect_that(res$Lb.Mean[1], equals(mean(c(862/3005,440/1530,
                                            506/3005,304/1530,
                                            914/3005,439/1530,
                                            723/3005,347/1530,
                                            1219/3363,638/1750,
                                            619/3363,309/1750,
                                            766/3363,402/1750,
                                            759/3363,401/1750,
                                            592/1760,
                                            743/1760,
                                            425/1760))))
  
  
  # Check result: Heterozygous balance standard deviation.
  expect_that(res$Lb.Stdv[1], equals(sd(c(862/3005,440/1530,
                                          506/3005,304/1530,
                                          914/3005,439/1530,
                                          723/3005,347/1530,
                                          1219/3363,638/1750,
                                          619/3363,309/1750,
                                          766/3363,402/1750,
                                          759/3363,401/1750,
                                          592/1760,
                                          743/1760,
                                          425/1760))))
  
  # Check result: 10 percentile for heterozygous balance.
  expect_that(res$Lb.Perc.10[1], equals(as.numeric(quantile(c(862/3005,440/1530,
                                                              506/3005,304/1530,
                                                              914/3005,439/1530,
                                                              723/3005,347/1530,
                                                              1219/3363,638/1750,
                                                              619/3363,309/1750,
                                                              766/3363,402/1750,
                                                              759/3363,401/1750,
                                                              592/1760,
                                                              743/1760,
                                                              425/1760), 0.10))))
  
  
})