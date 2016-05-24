context("Selection evaluation")

test_that("a single list of parameter specification is evaluated correctly", {
    climate <- data.frame(
      year = rep(1950:2009, each = 12),
      month = rep(1:12, 60),
      temp = rep(-10 * cos(seq(0, 2*pi, length.out = 12)), 60),
      prec = rep(seq(100, 220, length.out = 12), 60)
      )
    class(climate) <- c("tcclimate", "data.frame")
    chrono <- data.frame(rnorm(100))
    rownames(chrono) <- 1901:2000
    truncated_input <- truncate_input(chrono, climate, NULL, 1, FALSE)
    pmat <- make_pmat(truncated_input$climate)

    expect_that(eval_selection(pmat, 1:2)$month$single,
                equals(c("JAN", "FEB")))
    expect_that(eval_selection(pmat, 1:2)$aggregate[,1],
                equals(rep(-10, 50)))
    expect_that(round(eval_selection(pmat, .mean(-3:4))$aggregate[,1],
                      3),
                equals(rep(-0.909, 50)))
    expect_that(eval_selection(pmat, .sum(-11:3))$month$single,
                equals(c("nov", "dec", "JAN", "FEB", "MAR")))
})

test_that("the design matrix is constructed correctly", {
  climate <- data.frame(
    year = rep(1950:2009, each = 12),
    month = rep(1:12, 60),
    temp = rep(-10 * cos(seq(0, 2*pi, length.out = 12)), 60),
    prec = rep(seq(100, 220, length.out = 12), 60)
    )
  class(climate) <- c("tcclimate", "data.frame")
  chrono <- data.frame(rnorm(100))
  rownames(chrono) <- 1901:2000
  truncated_input <- truncate_input(chrono, climate, NULL, 1, FALSE)
  pmat <- make_pmat(truncated_input$climate)
  
  test_that(tc_design(.range(1:3), pmat)$aggregate[,4],
            equals(rep(-10, 50)))
  test_that(tc_design(.mean(1:3), pmat)$names[1],
            equals("prec.curr.jan"))
})
