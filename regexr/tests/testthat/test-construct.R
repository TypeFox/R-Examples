context("Checking construct")

test_that("construct makes a regular expression of class regexr with appropriate attributes",{

    m <- construct(
        space =
            "\\s+"
                %:)%"I see",
    
        simp =
            "(?<=(foo))",
    
        or =
            "(;|:)\\s*"
                %:)%"comment on what this does",
    
        "[a]s th[atey]"
    )

    m2 <- structure("\\s+(?<=(foo))(;|:)\\s*[a]s th[atey]", class = c("regexr", 
        "character"), subs = structure(list(space = structure("\\s+", class = c("subcom", 
        "character"), comment = "I see"), simp = "(?<=(foo))", or = structure("(;|:)\\s*", class = c("subcom", 
        "character"), comment = "comment on what this does"), "[a]s th[atey]"), .Names = c("space", 
        "simp", "or", "")), comments = structure(list(space = "I see", 
            simp = NULL, or = "comment on what this does", NULL), .Names = c("space", 
        "simp", "or", "")))

    expect_equivalent(m, m2)

})




