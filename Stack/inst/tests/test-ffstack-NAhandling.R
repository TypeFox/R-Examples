context("ffStack NA handling")


testdf1 <- data.frame(facNA.int=factor(c(letters[1:3],NA)),
                      intNA.fac=c(1,2,3,NA),
                      facNA.intNA=factor(c(letters[1:3],NA)),
                      fac.NA=factor(1),
                      allNA=NA,
                      emptystring=factor(c("",1:3)))

testdf2 <- data.frame(facNA.int=2,
                      intNA.fac=factor(letters[1:4]),
                      facNA.intNA=c(1,2,3,NA),
                      NA.fac=factor(2),
                      allNA=NA,
                      emptystring=factor(c("",4:6)))
ffdf1 <- ff::as.ffdf(.preparedf(testdf1))
ffdf2 <- ff::as.ffdf(.preparedf(testdf2))
stacked1 <- try(ffStack(ffdf2, testdf1))
stacked2 <- try(ffStack(ffdf1, testdf2))
stacked3 <- ffStack(ffdf1,ffdf2)
test_that("Stack() handles mixed factor/integers with NAs", {
    expect_true(nrow(stacked1)>0)
    expect_true(nrow(stacked2)>0)
})

testdf1 <- data.frame(facNA=(factor(c(1:4,NA))))
testdf2 <- data.frame(NAfac=factor(letters[1:5]),
                      NAlevel=addNA(factor(c(letters[23:26],NA))))
ffdf1 <- ff::as.ffdf(testdf1)
stacked <- ffStack(ffdf1,testdf2)
stacked <- as.data.frame(ffStack(ffdf1,testdf2))
nl <- vapply(stacked,function(x) length(levels(x)) , integer(1))
nna <- vapply(stacked,function(x) sum(is.na(x[])), integer(1))
test_that("Stack() handles mixed factor/integers with NAs", {
    expect_equal(nl, c(facNA=4,NAfac=5,NAlevel=4))
    expect_equal(nna, c(facNA=6,NAfac=5,NAlevel=6))
})


