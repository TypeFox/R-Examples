context("Checking comments")

test_that("comments grabs component chunks comments",{

    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
        "a", "b", structure("c", class = c("subcom", "character"), 
        comment = "A love note to your future self")), comments = list(
        NULL, NULL, "A love note to your future self"))

    
    x1 <- list(NULL, NULL, "A love note to your future self")

    expect_equivalent(comments(minimal), x1)

})

test_that("comments are changed by assigment",{
    
    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
        "a", "b", structure("c", class = c("subcom", "character"), 
        comment = "A love note to your future self")), comments = list(
        NULL, NULL, "A love note to your future self"))

    
    comments(minimal)[2] <- "\\s+[A-Z]|[0-9]"
    
    x2 <- list(NULL, "\\s+[A-Z]|[0-9]", "A love note to your future self")
    
    expect_equivalent(comments(minimal), x2)
})

test_that("comments are changed by setting (similar to `setNames`)",{
    
    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", "c"), comments = list(NULL, NULL, NULL))
    
    out <- set_comments(minimal, c("a", "b", "c"))
    
    out_check <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", "c"), comments = c("a", "b", "c"))
    
    
    expect_equivalent(out, out_check)

})