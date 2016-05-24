set.seed(8675309)

TestGenerateFactorData <- function() {
  foo <- GenerateFactorData(list(a = c("foo", "bar", "baz"),
                                 b = c("larry", "moe", "curly", "shemp")),
                            50)

  ## Check that we got a data frame back.
  checkTrue(is.data.frame(foo))

  ## It should have 50 rows.
  checkEquals(50, nrow(foo))

  ## All variables must be factore.
  checkTrue(all(sapply(foo, is.factor)))

  ## We should get all the levels, even if some of them don't appear
  ## because of random sampling.
  bar <- GenerateFactorData(list(a = c("foo", "bar", "baz"),
                                 b = c("larry", "moe", "curly", "shemp")),
                            2)
  checkEquals(4, nlevels(bar[, 2]))
  checkEquals(3, nlevels(bar[, 1]))
}
