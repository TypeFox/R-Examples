context("Checking get_construct")

test_that("get_construct creates a construct script that's character when `file = NULL`",{

    x <- structure("\\w*", class = c("regexr", "reverse_construct", "character"
        ), subs = structure(list(`1` = "\\w*"), .Names = "1"), comments = structure(list(
        `1` = "word characters (a-z, A-Z, 0-9, _) (0 or more times (matching the most amount possible))"), .Names = "1"), reverse_construct = structure(
        "construct(\n    `1` = \"\\\\w*\"               %:)%  \"word characters (a-z, A-Z, 0-9, _) (0 or more times (matching the most amount possible))\"\n)\n", class = c("reverse_construct", 
        "character"))
    )
      
    out <- "construct(\n    `1` = \"\\\\w*\"               %:)%  \"word characters (a-z, A-Z, 0-9, _) (0 or more times (matching the most amount possible))\"\n)\n"

    expect_equivalent(get_construct(x, file=NULL), out)
})

test_that("get_construct creates a construct script that's pretty printed when `file = \"\"`",{

    x <- structure("\\w*", class = c("regexr", "reverse_construct", "character"
        ), subs = structure(list(`1` = "\\w*"), .Names = "1"), comments = structure(list(
        `1` = "word characters (a-z, A-Z, 0-9, _) (0 or more times (matching the most amount possible))"), .Names = "1"), reverse_construct = structure(
        "construct(\n    `1` = \"\\\\w*\"               %:)%  \"word characters (a-z, A-Z, 0-9, _) (0 or more times (matching the most amount possible))\"\n)\n", class = c("reverse_construct", 
        "character"))
    )
    
    out2 <- c(
        "construct(", 
        "    `1` = \"\\\\w*\"               %:)%  \"word characters (a-z, A-Z, 0-9, _) (0 or more times (matching the most amount possible))\"", 
        ")"
    )


    expect_equivalent(capture.output(get_construct(x)), out2)
    
})