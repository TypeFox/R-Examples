context("Stack()")



testdf1 <- data.frame(both.int=1:4,
                      expand.factor=c("blue", "yellow"),
                      mixed.fac.int=factor(letters[1:4]),
                      date=as.Date("1983-11-22"),
                      df1only=rep(1:2, each=2),
                      mixed.fac.chr=I(c("a","b","NA",NA)))
testdf2 <- data.frame(both.int=5:24,
                      expand.factor=factor(rep(c(1:4, NA), 4)),
                      mixed.fac.int=1:4,
                      date=as.Date("1981-09-24"),
                      df2only=factor(c("c", "d")),
                      mixed.fac.chr=c("a","b",NA,"c"))

levels(testdf2$mixed.fac.chr) <- letters #overleveled
## put levels in a different order than sort(levels) would do, but
## don't make it an ordered factor. Result needs to be ordered
## to preserve this.
levels(testdf2$expand.factor) <- c("green", "blue", "red", "yellow")


test_that("Stack() warns", {
    expect_that(Stack(testdf2, testdf1), gives_warning())
    expect_that(Stack(testdf1, testdf2), gives_warning())
    expect_that(Stack(testdf1[,c(1,4)], testdf2[,c(1,4)]), does_not_give_warning())
})

test_that("Stack is quiet when it should be quiet", {
    expect_that(Stack(testdf2,testdf2), does_not_give_warning()) ## no longer warns for order
    expect_that(Stack(testdf1,testdf1), does_not_give_warning())
})

stacked1 <- try(Stack(testdf2, testdf1))
stacked2 <- try(Stack(testdf1, testdf2))
stacked1.chr <- try(Stack(testdf1,testdf2, mixed="c"))


test_that("Stack() stacks", {
    expect_true(class(stacked1)!="try-error")
    expect_true(class(stacked2)!="try-error")
    expect_equal(nrow(stacked1), nrow(testdf1) + nrow(testdf2))
    expect_equal(nrow(stacked2), nrow(testdf1) + nrow(testdf2))
})

test_that("Stack handles mixed factor/int/chr, new levels, and NAs", {
    expect_equal(levels(stacked1$mixed.fac.int), letters[1:4])
    expect_equal(sum(is.na(stacked1$df1only)), 20)
    expect_equal(sum(is.na(stacked1$df2only)), 4)
    expect_equal(sum(is.na(stacked1$expand.factor)), 4)
    expect_is(stacked1$mixed.fac.chr, "factor")
    expect_is(stacked1.chr$mixed.fac.chr, "character")

})

test_that("Factor level expansion works, with expected (non-alpha) order", {
    expect_equal(length(levels(stacked1$expand.factor)), 4)
    expect_equal(levels(stacked2$expand.factor),
                 c("green","blue","red","yellow"))
    expect_equal(length(levels(stacked1$expand.factor)),
                 length(unique(levels(stacked1$expand.factor))))
})

## test_that("Mixed factor-character with 'overleveled' factor drops levels",{
##     expect_
## }

test_that("Stack order doesn't matter", {
    expect_equivalent(stacked1[order(stacked1$both.int), names(stacked1)],
                      stacked2[order(stacked2$both.int), names(stacked1)])
})


testdf1 <- data.frame(facNA.int=factor(c(letters[1:3],NA)),
                      intNA.fac=c(1,2,3,NA),
                      facNA.intNA=factor(c(letters[1:3],NA)),
                      fac.NA=factor(1))
testdf2 <- data.frame(facNA.int=1,
                      intNA.fac=factor(letters[1:4]),
                      facNA.intNA=c(1,2,3,NA),
                      NA.fac=factor(2))
stacked1 <- try(Stack(testdf2, testdf1))
stacked2 <- try(Stack(testdf1, testdf2))
test_that("Stack() handles mixed factor/integers with NAs", {
    expect_true(nrow(stacked1)>0)
    expect_true(nrow(stacked2)>0)
})

testdf1 <- data.frame(facNA=(factor(c(1:4,NA))))
testdf2 <- data.frame(NAfac=factor(letters[1:5]),
                      NAlevel=addNA(factor(c(letters[23:26],NA))))
stacked <- Stack(testdf1,testdf2)
nl <- vapply(stacked,function(x) length(levels(x)) , integer(1))
nna <- vapply(stacked,function(x) sum(is.na(x)), integer(1))
test_that("Stack() handles mixed factor/integers with NAs", {
    expect_equal(nl, c(facNA=4,NAfac=5,NAlevel=4))
    expect_equal(nna, c(facNA=6,NAfac=5,NAlevel=6))
})


