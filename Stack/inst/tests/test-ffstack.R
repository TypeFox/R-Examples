context("ffStack()")

testdf1 <- data.frame(caseid=1:4,
                      expand.factor=factor(c("blue", "yellow")),
                      mixed.fac.int=factor(letters[1:4]),
                      date=as.Date("1983-11-22"),
                      df1only=rep(1:2, each=2),
                      mixed.fac.chr=c("a","b","NA",NA),
                      stringsAsFactors=FALSE)
ffdf1 <- ff::as.ffdf(.preparedf(testdf1))


testdf2 <- data.frame(caseid=5:24,
                      expand.factor=factor(rep(c(1:4, NA))),
                      mixed.fac.int=1:4,
                      date=as.Date("1981-09-24"),
                      df2only=factor(c("c", "d")),
                      mixed.fac.chr=c("a","b",NA,"c"))

levels(testdf2$mixed.fac.chr) <- letters #overleveled
## put levels in a different order than sort(levels) would do, but
## don't make it an ordered factor. Result needs to be ordered
## to preserve this.
levels(testdf2$expand.factor) <- c("green", "blue", "red", "yellow")
ffdf2 <- ff::as.ffdf(.preparedf(testdf2))


test_that("Stack() warns", {
    expect_that(ffStack(ffdf2, testdf1), gives_warning())
    expect_that(ffStack(ffdf1, testdf2), gives_warning())
    expect_that(ffStack(ffdf1[,c(1,4)], testdf2[,c(1,4)]), does_not_give_warning())
    #expect_that(ffStack(ffdf2,testdf2,verbose=TRUE), gives_warning()) ## warns for order (um, no?)
})
test_that("Stack is quiet when it should be quiet", {
    expect_that(ffStack(ffdf1,ffdf1), does_not_give_warning())
})

stacked1 <- try(ffStack(ffdf2, testdf1))
stacked2 <- try(ffStack(ffdf1, testdf2))
sdf1 <- as.data.frame(stacked1[order(stacked1$caseid[]), names(stacked1)])
sdf2 <- as.data.frame(stacked2[order(stacked2$caseid[]), names(stacked1)])

## Not applicable for ff stack
##stacked1.chr <- try(Stack(testdf1,testdf2, mixed="c"))

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
    expect_true(ff::is.factor(stacked1$mixed.fac.chr))

})

test_that("Factor level expansion works, with expected (non-alpha) order", {
    expect_equal(length(levels(stacked1$expand.factor)), 4)
    expect_equal(levels(stacked2$expand.factor),
                 c("green","blue","red","yellow"))
    expect_equal(length(levels(stacked1$expand.factor)),
                 length(unique(levels(stacked1$expand.factor))))
})

test_that("Stack order doesn't matter", {
    expect_equivalent(sdf1, sdf2)
})

## test upgrade from byte to short.
foo <- .compact(ffdf(foo=ff(c(-125,125))))
bar <- data.frame(foo=150)
baz <- ffStack(foo,bar)

foo <- .compact(ff::ffdf(foo=ff(as.Date(format(Sys.Date()))+1:5)))
bar <- ffStack(foo,data.frame(foo=foo$foo[]))

## TODO: ffStack is not overzealously expanding sizings
##       i.e., that it respects byte if byte is ok.
