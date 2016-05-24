test_that(
  "test.has_attributes.struct_with_attrs.returns_true_if_attr_exists",
  {
    x <- structure(list(a = 1), foo = 1, bar = 2)
    attrs <- c("names", "foo", "bar", "baz", NA)
    expected <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
    names(expected) <- attrs
    class(expected) <- c("vector_with_cause", "logical")
    # Note that base::attr treats missing attributes in the same way as
    # non-existent attributes (i.e., it returns NULL).  This should
    # behave the same way.
    assertive.base::cause(expected) <- c("", "", "", "no attr", "no attr")
    expect_equal(
      has_attributes(x, c("names", "foo", "bar", "baz", NA)),
      expected
    )
  }
)
