context("Boxplot")

expect_traces <- function(gg, n.traces, name) {
  stopifnot(is.ggplot(gg))
  stopifnot(is.numeric(n.traces))
  save_outputs(gg, paste0("boxplot-", name))
  L <- gg2list(gg)
  all.traces <- L$data
  no.data <- sapply(all.traces, function(tr) {
    is.null(tr[["x"]]) && is.null(tr[["y"]])
  })
  has.data <- all.traces[!no.data]
  expect_equal(length(has.data), n.traces)
  list(traces=has.data, layout=L$layout)
}

test_that("geom_boxplot gives a boxplot", {
  gg <- ggplot(mtcars, aes(factor(cyl), mpg)) + geom_boxplot()
  
  L <- save_outputs(gg, "boxplot")

  # right nb. traces
  expect_equal(length(L$data), 1)
  # right type for 1st trace
  expect_identical(L$data[[1]]$type, "box")
})

test_that("you can make a boxplot for a distribution of datetimes", {
  dist <- c(10, 20, 33, 40, 11, 12, 11)
  dist <- as.POSIXct(paste0("2014-09-19 10:00:", dist))
  df <- data.frame(y=dist)
  df$x <- 0
  
  bp <- ggplot(df) + geom_boxplot(aes(x, y))
  
  L <- save_outputs(bp, "boxplot-datetime")
  
  expect_equal(length(L$data), 1)  # 1 trace
  expect_identical(L$data[[1]]$type, "box")
  expect_identical(L$data[[1]]$y, as.numeric(df$y))
})

# check legend shows up when each box-and-whiskers has a fill
# make ggplot2
m <- ggplot(mtcars, aes(factor(cyl), mpg))
p <- m + geom_boxplot(aes(fill = factor(cyl)))
# tests
test_that("legends for boxplot", {
  info <- expect_traces(p, 3, "legends_for_fill")
  tr <- info$traces
  la <- info$layout
  expect_identical(tr[[1]]$type, "box")
  # check legend exists
  expect_identical(la$showlegend, TRUE)
  # check legend for each fill exists
  for (i in 1:3) {
    expect_identical(tr[[i]]$showlegend, TRUE)
  }
  # check the fill colors are correct
  g <- ggplot_build(p)
  fill.colors <- unique(g$data[[1]]["fill"])[,1]
  for (i in 1:3) {
    plotly.color <- as.integer(strsplit(gsub("[\\(\\)]|rgb", "", 
                tr[[i]]$fillcolor), split = ",")[[1]])
    names(plotly.color) <- c("red", "green", "blue")
    expect_equal(plotly.color, col2rgb(fill.colors[i])[,1], 
                 tolerance = 1)
  }
})

dat <- data.frame(
  cond = factor(rep(c("A", "B", "C", "D"), each = 200)), 
  col = factor(rep(c("C1", "C2"), each = 400)), 
  rating = c(rnorm(200), rnorm(200, mean=.8), rnorm(200, mean=.4), rnorm(200, mean=.2))
)
g <- ggplot(dat, aes(x = cond, y = rating)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = col))

test_that("correct # of unique fillcolors", {
  L <- save_outputs(g, "boxplot-fillcolor")
  expect_equal(length(L$data), 2)
  expect_identical(L$data[[1]]$type, "box")
  fills <- sapply(L$data, "[[", "fillcolor")
  expect_equal(length(unique(fills)), length(unique(dat$col)))
})
