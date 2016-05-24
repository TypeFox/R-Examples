test_that(
  "All functions from sub packages are exported",
  {
    pkgs <- paste0(
      "assertive.", 
      c(
        "base", "properties", "types",
        "numbers", "strings", "sets",
        "matrices", "models", "datetimes",
        "data", "data.uk", "data.us",
        "reflection", "code", "files"
      )
    )
    
    x <-
      unlist(
        lapply(
          pkgs,
          function(pkg)
          {
            library(pkg, character.only = TRUE)
            ls(paste0("package:", pkg))
          }
        )
      )
      library(assertive)
      y <- ls("package:assertive")
      
      should_be_exported_but_isnt <- setdiff(x, y)
      shouldnt_exist <- setdiff(y, x)
      expect_equal(
        length(should_be_exported_but_isnt), 
        0, 
        label = toString(should_be_exported_but_isnt)
      )
      expect_equal(
        length(shouldnt_exist), 
        0,
        label = toString(shouldnt_exist)
      )
  }
)

