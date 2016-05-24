context("YahooInPlaySource")

test_that("YahooInPlaySource",{
	
	minlengthcorp <- 1
		
	testcorp <- WebCorpus(YahooInplaySource())
	lengthcorp <- length(testcorp)
	# Check Corpus object
	expect_that(length(testcorp) >= minlengthcorp, is_true())
	expect_that(class(testcorp), equals(c("WebCorpus","VCorpus","Corpus")))
	
	# Check Content
	#expect_that(all(sapply(testcorp, nchar) > 0), is_true())
	contentlength <- sapply(testcorp, function(x) 
				if( length(content(x)) < 1) 0 else nchar(content(x)))	
	contentratio <- length(which(contentlength > 0)) / length(testcorp)
	expect_that(contentratio > 0.5, is_true())
	
	# Check Meta Data
	datetimestamp <- lapply(testcorp, function(x) meta(x, "datetimestamp"))
	#FIXME: Date should be fixed
	expect_that(all(sapply(datetimestamp, function(x) class(x)[1] == "character")), is_true())
	
	heading <- lapply(testcorp, function(x) meta(x, "heading")[1])
	expect_that(all(sapply(heading, function(x) class(x)[1] == "character")), is_true())
	expect_that(all(sapply(heading, nchar) > 0), is_true())
	
	id <- lapply(testcorp, function(x) meta(x, "id")[1])
	expect_that(all(sapply(id, function(x) class(x)[1] == "character")), is_true())
	expect_that(all(sapply(id, nchar) > 0), is_true())
	
	testcorp <- testcorp[1:length(minlengthcorp)]
	testcorp <- corpus.update(testcorp)
	expect_that(length(testcorp) >= lengthcorp, is_true())
	
	cat(" | Contentratio: ", sprintf("%.0f%%", contentratio * 100))
})

