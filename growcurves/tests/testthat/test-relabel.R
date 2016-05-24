## library(growcurves)


##
## perform tests of relabel function
##

context("Relabel Function Accuracy")

##
## Load simulation dataset without nuisance covariates
## (Two treatment levels, {0,1}, and no nuisance covariates)
data(datsim)

test_that("numeric label input correctly handled", {

	lab 	<- relabel(label.input = datsim$trt, start = 0)

	expect_that(length(names(lab)), equals(3)) ## returns correct number of objects
	expect_that(lab$label.new, is_a("integer")) ## test that is.integer(label.new) == TRUE
	expect_that(which(lab$labeli.u == min(lab$labeli.u)), equals(1)) ## sort order of labeli.u
	expect_that(sum(lab$label.new - datsim$trt),is_equivalent_to(0)) ## check order of label.new same as input label
})

test_that("character label input correctly handled", {
	
	trt.lab <- factor(datsim$trt,labels=c("fred","barney"))
	lab  <- relabel(label.input = trt.lab, start = 0)

	expect_that(factor(lab$labeli.u), is_equivalent_to(unique(trt.lab))) ## test that labeli.u captures unique input values
	expect_that(length(unique(lab$label.new)), is_equivalent_to(length(lab$labeli.u))) ## test label.new has right number of unique values
      	expect_that(lab$label.new, is_a("integer")) ## test that is.integer(label.new) == TRUE
	expect_that(min(lab$label.new), equals(0)) ## test that minimum value in label.new is 0

})