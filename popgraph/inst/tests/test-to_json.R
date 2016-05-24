context("as_json.R")

test_that("tests", {
  a <- matrix( 0,nrow=4,ncol=4)
  a[1,2] <- a[1,3] <- a[2,3] <- a[1,4] <-1
  a <- a + t(a)
  graph <- as.popgraph(a)
  
  json <- to_json(graph)  
  expect_that( json, is_a("character"))
  expect_that( json, is_equivalent_to("var myjson = '{ \"nodes\":[{\"name\":\"node-1\",\"group\":\"All\"}, {\"name\":\"node-2\",\"group\":\"All\"}, {\"name\":\"node-3\",\"group\":\"All\"}, {\"name\":\"node-4\",\"group\":\"All\"}], \"links\":[{\"source\":0,\"target\":1}, {\"source\":0,\"target\":2}, {\"source\":0,\"target\":3}, {\"source\":1,\"target\":2}]}';"))
  
  expect_that( to_json(FALSE), throws_error() )
  
})