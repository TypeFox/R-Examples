if(require(testthat)) {
  require(traitr)

  context("model")
  i <- Model$proto()

  ## init
  i$prop1 <- "property1"
  i$init()
  expect_that(i$prop1 == "property1", is_true())
  expect_that(i$get_prop1() == "property1", is_true())
  i$set_prop1("new prop")
  expect_that(i$get_prop1() == "new prop", is_true())

  
  ## observers
  i$prop2 <- "property 2"; i$init()
  i$property_prop1_value_changed <- function(.,value, old_value, ...) .$set_prop2(value)
  i$set_prop1("new value")
  expect_that(i$get_prop2() == "property 2", is_true())
  i$add_observer(i)
  i$set_prop1("another value")
  expect_that(i$get_prop2() == "another value", is_true())

  ## remove observer
  i$remove_observer(i)
  i$set_prop1("yet another value")
  expect_that(i$get_prop2() == "yet another value", is_false())

## some more -- trim down
    m <- aModel(a=1, b= 2)

  ## get
  expect_that(m$get_a() == 1, is_true())

  ## setter/getter
  m$set_a(2)
  expect_that(m$get_a() == 2, is_true())

  ## observe self
  m$property_a_value_changed <- function(., ...) .$set_b(.$get_a())
  m$set_a(3)
  expect_that(m$get_b() == 3, is_false())
  m$add_observer(m)

  m$set_a(4)
  expect_that(m$get_b() == 4, is_true())


}
