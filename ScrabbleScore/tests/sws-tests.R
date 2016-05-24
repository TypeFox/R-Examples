library(testthat)
library(ScrabbleScore)

context("sws")

test_that("score correct for single words sws where distribution doesn't matter",{
  expect_equal(sws("the"),6)
  expect_equal(sws("quick"),20)
  expect_equal(sws("brown"),10)
  expect_equal(sws("fox"),13)
  expect_equal(sws("jumps"),16)
  expect_equal(sws("over"), 7) 
  expect_equal(sws("lazy"), 16)
  expect_equal(sws("dog"),5)
})

test_that("score correct for vectors of words on words where distribution doesn't matter",{
  expect_equal(sws(c("the","quick","brown","fox","jumps","over","lazy","dog")),c(6,20,10,13,16,7,16,5))
})

test_that("score correct for single words with only 4 tiles",{
  expect_equal(sws("ddddddun"),10) 
  expect_equal(sws("huuuuuuuuuurr"),10)
  expect_equal(sws("ssssssssss"),4)
  expect_equal(sws("llllllollllll"),5)
})

test_that("score correct for vector of words with only 4 tiles",{
  expect_equal(sws(c("ddddddun","huuuuuuuuuurr","ssssssssss","llllllollllll")),c(10,10,4,5))
})

test_that("score correct for single words with only 3 tiles",{
  expect_equal(sws("gggggggrrr") ,(3*sls("g"))+(3*sls("r"))) 
})

test_that("score correct for single words with only 2 tiles",{
  expect_equal(sws("bbbbb"),(2*sls("b")))
  expect_equal(sws("cccp"),(2*sls("c")+sls("p")))   
  expect_equal(sws("mmmmmmmm"),(2*sls("m")))
  expect_equal(sws("ppppppplease"),(2*sls("p")+sls("l")+2*sls("e")+sls("a")+sls("s"))) 
  expect_equal(sws("fffffudge"), 2*sls("f")+sls("u")+sls("d")+sls("g")+sls("e"))
  expect_equal(sws("hhhhhhhmmmmmm") , 2*sls("m") + 2*sls("h") )      
  expect_equal(sws("vvvvvvvvvv") , 2*sls("v") )
  expect_equal(sws("wwwwwwwhat"), 2*sls("w")+sls("h")+sls("a")+sls("t"))
  expect_equal(sws("whyyyyyyyyy") , 2*sls("y")+sls("w")+sls("h"))
})

test_that("score correct for vector of words with only 2 tiles",{
  expect_equal(sws(c("bbbbb","cccp","mmmmmmmm","ppppppplease","fffffudge","hhhhhhhmmmmmm","vvvvvvvvvv","wwwwwwwhat","whyyyyyyyyy")),
                   c((2*sls("b")),(2*sls("c")+sls("p")),(2*sls("m")),(2*sls("p")+sls("l")+2*sls("e")+sls("a")+sls("s")),
                     2*sls("f")+sls("u")+sls("d")+sls("g")+sls("e"), 2*sls("m") + 2*sls("h"), 2*sls("v"), 2*sls("w")+sls("h")+sls("a")+sls("t"),
                     2*sls("y")+sls("w")+sls("h")))
})


test_that("score correct for single words with only 1 tile",{
  expect_equal(sws("xxx"),sls("x"))
  expect_equal(sws("kicked"),sls("k")+sls("i")+sls("c")+sls("e")+sls("d"))
  expect_equal(sws("jjzzz"),sls("z")+sls("j"))
  expect_equal(sws("quiq"),sls("q")+sls("u")+sls("i"))
  expect_equal(sws("zzz"),sls("z"))
})

test_that("score correct for vector of words with only 2 tiles",{
  expect_equal(sws(c("xxx","kicked","jjzzz","quiq","zzz")),
               c(sls("x"),sls("k")+sls("i")+sls("c")+sls("e")+sls("d"),sls("z")+sls("j"),sls("q")+sls("u")+sls("i"),sls("z")))
})

test_that("the parameter 'check.valid' if TRUE returns score if possilbe otherwise 0",{
  expect_equal(sws("zzz",check.valid=TRUE),10)
  expect_equal(sws("zzzz",check.valid=TRUE),0)
  expect_equal(sws(c("zzz","zzzz"),check.valid=TRUE),c(10,0))
})

test_that("the parameters 'only.possible' if FALSE will provide the correct impossible scores",{
  expect_equal(sws("zzz",only.possible=FALSE),30)
  expect_equal(sws(c("zzz","zzzz"),only.possible=FALSE),c(30,40))
  
})
