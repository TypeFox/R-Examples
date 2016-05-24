context("Checking %:)%")

test_that("%:)% adds comments",{

    minimal <- construct("a", "b", "c" %:)% "A love note to your future self")
    
    x1 <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", structure("c", class = c("subcom", "character"), 
                comment = "A love note to your future self")), comments = list(
            NULL, NULL, "A love note to your future self"))
    
    expect_equivalent(x1, minimal)
    
    x2 <- structure(list(class = c("regexr", "character"), subs = list(
            "a", "b", structure("c", class = c("subcom", "character"), 
                comment = "A love note to your future self")), 
            comments = list(NULL, NULL, "A love note to your future self")), 
        .Names = c("class", 
        "subs", "comments"))
    
    
    expect_equivalent(x2, attributes(minimal))
    
    x3 <- list(NULL, NULL, "A love note to your future self")

    expect_equivalent(x3, comments(minimal))

})

test_that("%:)% & %comment% are equal",{

    expect_equivalent("A" %:)% "B", "A" %comment% "B")
})

