context("Checking subs")

test_that("subs grabs component chunks",{

    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", structure("c", class = c("subcom", "character"), 
                comment = "A love note to your future self")), comments = list(
            NULL, NULL, "A love note to your future self"))
    
    x1 <- list("a", "b", structure("c", class = c("subcom", "character"), 
        comment = "A love note to your future self"))
    
    expect_equivalent(subs(minimal), x1)

})

test_that("subs is changed by assigment",{
    
    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", structure("c", class = c("subcom", "character"), 
                comment = "A love note to your future self")), comments = list(
            NULL, NULL, "A love note to your future self"))
    
    subs(minimal)[2] <- "\\s+[A-Z]|[0-9]"
    
    x2 <- list("a", "\\s+[A-Z]|[0-9]", structure("c", class = c("subcom", 
        "character"), comment = "A love note to your future self"))
    
    expect_equivalent(subs(minimal), x2)
})

test_that("subs are changed by setting (similar to `setNames`)",{
    
    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", "c"), comments = list(NULL, NULL, NULL))
    
    out <- set_subs(minimal, c("\\s+", "(?:foo)", "[A-Za-z-]*"))
    
    out_check <- structure("\\s+(?:foo)[A-Za-z-]*", class = c("regexr", "character"
        ), subs = c("\\s+", "(?:foo)", "[A-Za-z-]*"), comments = list(
            NULL, NULL, NULL))
    
    expect_equivalent(out, out_check)

})

