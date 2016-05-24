context("Compare results with those in Baggelaar et al (2010)")

test_that("EQR-results on page 22 can be reproduced.", {
    data(EQR)
    result <- mya(EQR)
    expect_that(result$MYA, equals(0.167985609))
    expect_that(result$q05, equals(0.075350782))
    expect_that(result$q95, equals(0.333436559))
    expect_that(result$PROB_LTT, equals(0.988412826))
    expect_that(result$PROB_GTT, equals(0.011587173501))
})

test_that("DCA-results on page 24 can be reproduced.", {
    data(DCA)
    result <- mya(DCA)
    expect_that(result$MYA, equals(9.5))
    expect_that(result$q05, equals(7.5852374065))
    expect_that(result$q95, equals(11.4147625935))
    expect_that(result$PROB_LTT, equals(0.7372895))
    expect_that(result$PROB_GTT, equals(0.2627105))
})
