context("retornaValor return")

test_that("retornaValor o valor retornado", {
  expect_equal(retornaValor("Abc"), "Abc")
  expect_equal(retornaValor("02/2014"),"02/2014")
  expect_equal(retornaValor(01), 1)
})
