library(rjson);

context("Abstract Syntax Tree");

test_that("Basics",{
 expect_identical(
  toJSON(r()$db("a")$table("b")$add()$query),
  "[24,[[15,[[14,[\"a\"]],\"b\"]]]]");
 expect_identical(
  toJSON(r()$db("a")$table("b")$query),
  "[15,[[14,[\"a\"]],\"b\"]]");

});

test_that("Functions",{
 Q<-r()$funcall(function(x) x,777)$query;
 expect_identical(toJSON(Q),
  "[64,[[69,[[2,[1]],[10,[1]]]],777]]")
 expect_error(r()$filter(function() r()$add(1))$query);
 expect_error(r()$filter(list(a=function(x) x))$query,
  "Functions can only exist as direct term arguments.");
})

test_that("Make array appears",{
 Q<-r()$db(c("a","b","c"))$query;
 expect_identical(toJSON(Q),"[14,[[2,[\"a\",\"b\",\"c\"]]]]");
})

test_that("Expressions",{
 Q<-r()$add(r()$add(1,2),4)$query;
 expect_identical(toJSON(Q),"[24,[[24,[1,2]],4]]");
})
