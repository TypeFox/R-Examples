context("Checking unglue")

test_that("unglue extracts elements chunks from a regecr object",{

    minimal <- structure("abc", class = c("regexr", "character"), subs = list(
            "a", "b", structure("c", class = c("subcom", "character"), 
                comment = "A love note to your future self")), comments = list(
            NULL, NULL, "A love note to your future self"))
    
    x <- structure(list("a", "b", structure("c", class = c("subcom", "character"
        ), comment = "A love note to your future self")), class = c("unglued", 
        "list"))
    
    
    expect_equivalent(unglue(minimal), x)

})

test_that("unglue prints as expected",{

    minimal_unglue <- structure(list("a", "b", structure("c", class = c("subcom", "character"
        ), comment = "A love note to your future self")), class = c("unglued", 
        "list"))

    expect_equivalent(capture.output(regexr:::print.unglued(minimal_unglue)),
        c("[[1]]", "[1] \"a\"", "", "[[2]]", "[1] \"b\"", "", "[[3]]", 
        "[1] \"c\"", ""))

})
