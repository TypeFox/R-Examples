context("Checking F-test significance helper functions")

test_that("Can judge when an F-statistic is significant", {
    model.summary <- summary(lm(weight ~ feed, data = chickwts)) # F > 15
    expect_equal(IsFSignificant(model.summary), TRUE)
  }
)

test_that("Can judge when an F-statistic is NOT significant", {
    set.seed(1001)
    value <- rnorm(n = 100)
    group <- rep(1:4, times = 25)
    fake.data <- data.frame(value, group)
    model.summary <- summary(lm(value ~ group, data = fake.data))

    expect_equal(IsFSignificant(model.summary), FALSE)
  }
)

test_that("An F-statistic can be > 1 but still NOT significant", {
    set.seed(1001)
    red <- rnorm(n = 30)
    set.seed(1001)
    green <- rnorm(n = 30, mean = 0.3)
    set.seed(1001)
    blue <- rnorm(n = 30, mean = -0.3)

    value <- c(red, green, blue)
    group <- rep(c("red", "green", "blue"), each = 30)
    fake.data <- data.frame(group, value)
    model.summary <- summary(lm(value ~ group, data = fake.data)) # Yields an F > 1, but is not significant

    expect_equal(IsFSignificant(model.summary), FALSE)
  }
)