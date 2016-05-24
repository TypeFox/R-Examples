context("Checking names")

test_that("names are changed by assigment",{
    
    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
        "a", "b", structure("c", class = c("subcom", "character"), 
            comment = "A love note to your future self")), comments = list(
        NULL, NULL, "A love note to your future self"))

    
    names(minimal)[2:3] <- c("Foo", "Bar")
    
    x <- c("", "Foo", "Bar")
    
    expect_equivalent(names(minimal), x)
})


test_that("names are changed by setting (similar to `setNames`)",{
    
    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", "c"), comments = list(NULL, NULL, NULL))
    
    out <- set_names(minimal, 1:3)
    
    out_check <- structure("abc", class = c("regexr", "character"), subs = structure(list(
            `1` = "a", `2` = "b", `3` = "c"), .Names = c("1", "2", "3"
        )), comments = structure(list(`1` = NULL, `2` = NULL, `3` = NULL), .Names = c("1", 
        "2", "3")))
    
    expect_equivalent(out, out_check)
})
