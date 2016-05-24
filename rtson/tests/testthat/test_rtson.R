library(rtson)

context("TSON test")

expect_equal_list = function(found, expected){
  
  if ( is.null(found) != is.null(expected) ) stop(paste0("expect_equal_list length failed : ", label, " : expected ", expected , " found " , found))
  expect_equal(length(found), length(expected))
  if (length(found) > 0) {
    
    found_names = names(found)
    expected_names = names(expected)
    expect_equal(found_names, expected_names)
    
    if (length(found) > 1){
      for (i in seq(1,length(found))){
        expect_equal_list(found[[i]], expected[[i]])
      }
    } else {
      expect_equal(as.vector(found[[1]]), as.vector(expected[[1]]) )
    }
  }
}

test_that("Empty list", {
  list = list()
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Empty map", {
  list = list()
  bytes = toTSON(tson.map(list))
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Null scalar character", {
  list = list(a=tson.character(NULL))
  list = list(a=NULL)
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Simple list", {
  list = list(tson.character("a"), TRUE, FALSE, tson.int(42L), tson.double(42.0) )
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Simple int32 list", {
  list = as.integer(c(42,42))
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Simple cstring list", {
  list = c("42.0","42")
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Simple map", {
  list = list(a=tson.character("a"), i=tson.int(42L), d=tson.double(42.0) )
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("Simple map of int32, float32 and float64 list", {
  list = list(i=42L, f=tson.float32.vec(42.0) , d=as.double(42.0) )
  bytes = toTSON(list)
  #   print(as.integer(bytes))
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

test_that("All types", {
  list = list(string=tson.scalar("string"),
              integer=tson.scalar(42L),  
              double=tson.scalar(42),  
              bool=tson.scalar(TRUE),  
              cstringlist=c("42",42.0),
              uint8=tson.uint8.vec(c(42,0)),
              uint16=tson.uint16.vec(c(42,0)),
              uint32=tson.uint32.vec(c(42,0)),
              int8=tson.int8.vec(c(42,0)),
              int16=tson.int16.vec(c(42,0)),
              int32=as.integer(c(42,0)),
              float32=tson.float32.vec(c(0.0, 42.0)),
              float64=c(42.0,42.0),
              map=list(x=42, y=42, label=tson.scalar("mylabel")),
              list=list("42",42)
  )
  bytes = toTSON(list)
  object = fromTSON(bytes)
  expect_equal_list(object, list)
})

