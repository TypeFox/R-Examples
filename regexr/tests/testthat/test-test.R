context("Checking test")

test_that("test matches epected regular expression validity (is valid)",{

    m <- structure("\\s+(?<=(foo))(;|:)\\s*[a]s th[atey]", class = c("regexr", 

        "character"), subs = structure(list(space = structure("\\s+", class = c("subcom", 
        "character"), comment = "I see"), simp = "(?<=(foo))", or = structure("(;|:)\\s*", class = c("subcom", 
        "character"), comment = "comment on what this does"), "[a]s th[atey]"), .Names = c("space", 
        "simp", "or", "")), comments = structure(list(space = "I see", 
            simp = NULL, or = "comment on what this does", NULL), .Names = c("space", 
        "simp", "or", "")))
    
    x1 <- structure(list(subs = TRUE, subexpressions = structure(c(TRUE, TRUE, 
        TRUE, TRUE), .Names = c("space", "simp", "or", ""))), .Names = c("regex", 
        "subexpressions"))
    
    expect_equivalent(test(m), x1)
})

test_that("test matches epected regular expression validity (not valid)",{
    
    m2 <- structure("\\s+(?<=(foo))(;|:)\\s*[a]s th[atey](([A-Z]|(\\d{5}))", class = c("regexr", 
        "character"), subs = structure(list(space = structure("\\s+", class = c("subcom", 
        "character"), comment = "I see"), simp = "(?<=(foo))", or = structure("(;|:)\\s*", class = c("subcom", 
        "character"), comment = "comment on what this does"), "[a]s th[atey]", 
            "(", "([A-Z]|(\\d{5})", ")"), .Names = c("space", "simp", 
        "or", "", "", "", "")), comments = structure(list(space = "I see", 
            simp = NULL, or = "comment on what this does", NULL, NULL, 
            NULL, NULL), .Names = c("space", "simp", "or", "", "", "", 
        "")))
    
    expect_warning(test(m2))
    
    x2 <- structure(list(subs = FALSE, subexpressions = structure(c(TRUE, TRUE, 
        TRUE, TRUE, FALSE, FALSE, FALSE), .Names = c("space", "simp", 
        "or", "", "", "", ""))), .Names = c("regex", "subexpressions"))

    out <- test(m2, quiet=TRUE)
    expect_equivalent(out, x2)
})

