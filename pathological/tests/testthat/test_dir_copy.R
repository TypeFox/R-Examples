dir2 <- function(path = ".", pattern = NULL, recursive = TRUE, ...)
{
  dir(  
    path,
    pattern      = pattern, 
    recursive    = recursive, 
    all.files    = TRUE,
    full.names   = FALSE,
    include.dirs = TRUE,
    ...
  )
}

test_that(
  "copy_dir works with recursive = FALSE",
  {
    source_dir <- R.home("etc")
    target_dir <- file.path(tempdir(), "etc")
    on.exit(unlink(target_dir))
    copy_dir(source_dir, target_dir, recursive = FALSE)
    expect_equal(
      dir2(source_dir, recursive = FALSE), 
      dir2(target_dir, recursive = FALSE)
    )   
  }
)

test_that(
  "copy_dir works with recursive = TRUE",
  {
    source_dir <- R.home("etc")
    target_dir <- file.path(tempdir(), "etc")
    on.exit(unlink(target_dir, recursive = TRUE))
    copy_dir(source_dir, target_dir)
    expect_equal(dir2(source_dir), dir2(target_dir))   
  }
)

