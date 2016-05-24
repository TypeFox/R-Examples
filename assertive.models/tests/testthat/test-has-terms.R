test_that(
  "test.has_terms.without_terms.returns_false",
  {
    x <- 1:10
    actual <- has_terms(x)
    expect_false(actual)
    expect_equal(
      cause(actual), 
      noquote("x has no terms component nor attribute.")
    )
  }
)

test_that(
  "test.has_terms.with_terms_component.returns_false",
  {
    x <- list(terms = 1:10)
    expect_true(has_terms(x))
  }
)

test_that(
  "test.has_terms.with_terms_attribute.returns_false",
  {
    x <- 1:10
    attr(x, "terms") <- 1:10
    expect_true(has_terms(x))
  }
)

test_that(
  "test.has_terms.lm_model.returns_false",
  {
    x <- lm(uptake ~ conc, CO2)
    expect_true(has_terms(x))
  }
)
