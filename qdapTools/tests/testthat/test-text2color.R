context("Checking text2color")

test_that("text2color produce proper vectors and types",{

	    x <- structure(list(X1 = structure(c(3L, 1L, 8L, 4L, 7L, 2L, 2L, 2L,
        4L, 8L, 4L, 3L, 5L, 3L, 1L, 8L, 7L, 2L, 1L, 6L), .Label = c("a",
        "and", "in", "is", "of", "that", "the", "to"), class = "factor")),
        .Names = "X1", row.names = c(NA, -20L), class = "data.frame")
    
    expected1 <- c("white", "white", "white", "blue", "red", "green", "green", 
        "green", "blue", "white", "blue", "white", "white", "white", 
        "white", "white", "red", "green", "white", "white")
    
    expected2 <- structure(list(X1 = structure(c(3L, 1L, 8L, 4L, 7L, 2L, 2L, 2L, 
        4L, 8L, 4L, 3L, 5L, 3L, 1L, 8L, 7L, 2L, 1L, 6L), .Label = c("a", 
        "and", "in", "is", "of", "that", "the", "to"), class = "factor"), 
            X2 = c("white", "white", "white", "red", "red", "red", "red", 
            "red", "red", "white", "red", "white", "white", "white", 
            "white", "white", "red", "red", "white", "green")), .Names = c("X1", 
        "X2"), row.names = c(NA, -20L), class = "data.frame")
    
    
    
    expect_warning(text2color(x$X1, c("the", "and", "is"), c("red", "green", "blue")))
    expect_equivalent(
        text2color(x$X1, c("the", "and", "is"), c("red", "green", "blue", "white")),
        expected1
    )
    
    x$X2 <- text2color(x$X1, list(c("the", "and", "is"), "that"),
        c("red", "green", "white"))
    expect_equivalent(x, expected2)
    
	expect_true(is.character(text2color(x$X1, c("the", "and", "is"), 
		c("red", "green", "blue", "white"))))

})
