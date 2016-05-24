context("ml.arules")

test_that("ml.arules works", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
  mlBaskets <- ml.load.sample.data(myConn, "baskets", "baskets-test")
  db <- "rfml"
  expect_message(ml.add.index(x = mlBaskets$lineItem1productName, scalarType = "string", database =  db, conn = myConn), "Range element index created on productName")
  # We need to wait so that the index gets updated before using a function that leverage it
  Sys.sleep(10)
  itemsets <- ml.arules(mlBaskets, mlBaskets$lineItem1productName, support = 0.22, confidence = 0.01, target = "frequent itemsets")
  expect_is(itemsets, "itemsets")
  expect_equal(length(itemsets), 13)
  rules <- ml.arules(mlBaskets, mlBaskets$lineItem1productName, support = 0.22, confidence = 0.01)
  expect_is(rules, "rules")
  expect_equal(length(rules), 23)
  rm.ml.data.frame(mlBaskets)
})


