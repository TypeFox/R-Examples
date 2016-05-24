create_expected_decomposed_path <- function(dirname, filename, extension, row.names)
{
  structure(
    data.frame(
      dirname          = dirname,
      filename         = filename,
      extension        = extension,
      row.names        = row.names,
      stringsAsFactors = FALSE
    ), 
    class = c("decomposed_path", "data.frame")
  )
}

test_that(
  "decompose_path works with a zero length input",
  {
    x <- character()
    x2 <- NULL
    expected <- create_expected_decomposed_path(
      dirname   = character(), 
      filename  = character(), 
      extension = character(),
      row.names   = character()
    )
    expect_equal(decompose_path(x), expected)
    expect_equal(decompose_path(x2), expected)
  }
)

test_that(
  "decompose_path handles paths with no directory and a single extension in the filename.",
  {
    x <- "foo.tgz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = pwd,
      filename         = "foo",
      extension        = "tgz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles paths with a directory and a single extension in the filename.",
  {
    x <- "somedir/foo.tgz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = file.path(pwd, "somedir"),
      filename         = "foo",
      extension        = "tgz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles paths with no directory and a double extension in the filename.",
  {
    x <- "foo.tar.gz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = pwd,
      filename         = "foo",
      extension        = "tar.gz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles paths with a directory and a double extension in the filename.",
  {
    x <- "somedir/foo.tar.gz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = file.path(pwd, "somedir"),
      filename         = "foo",
      extension        = "tar.gz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)


test_that(
  "decompose_path handles paths with no directory and no extension in the filename.",
  {
    x <- "foo"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = pwd,
      filename         = "foo",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles paths with a directory and no extension in the filename.",
  {
    x <- "somedir/foo"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = file.path(pwd, "somedir"),
      filename         = "foo",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles filenames containing a '.' and an extension.",
  {
    x <- "foo. bar.zip"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = pwd,
      filename         = "foo. bar",
      extension        = "zip",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles backslashes in the directory name.",
  {
    x <- "somedir\\foo.tgz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = file.path(pwd, "somedir"),
      filename         = "foo",
      extension        = "tgz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles mixed forward and backslashes in the directory name.",
  {
    x <- "somedir\\another dir/foo.tgz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = file.path(pwd, "somedir", "another dir"),
      filename         = "foo",
      extension        = "tgz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles absolute paths to directories.",
  {
    x <- R.home()
    expected <- create_expected_decomposed_path(
      dirname          = normalizePath(R.home(), "/", mustWork = FALSE),
      filename         = "",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles '~'.",
  {
    x <- "~"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = normalizePath("~", "/", mustWork = FALSE),
      filename         = "",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles files inside '~'.",
  {
    x <- "~/foo.tgz"
    expected <- create_expected_decomposed_path(
      dirname          = normalizePath(dirname(x), "/", mustWork = FALSE),
      filename         = "foo",
      extension        = "tgz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles the current directory as '.'.",
  {
    x <- "."
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = pwd,
      filename         = "",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles the parent directory as '..'.",
  {
    x <- ".."
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = dirname(pwd),
      filename         = "",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles files inside '.'.",
  {
    x <- "./foo.tgz"
    pwd <- getwd()
    expected <- create_expected_decomposed_path(
      dirname          = pwd,
      filename         = "foo",
      extension        = "tgz",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles empty strings.",
  {
    x <- ""
    expected <- create_expected_decomposed_path(
      dirname          = "",
      filename         = "",
      extension        = "",
      row.names        = x
    )
    expect_equal(decompose_path(x), expected)
  }
)

test_that(
  "decompose_path handles missing paths.",
  {
    x <- NA
    expected <- create_expected_decomposed_path(
      dirname          = NA_character_,
      filename         = NA_character_,
      extension        = NA_character_,
      row.names        = "<NA>"
    )
    expect_warning(
      actual <- decompose_path(x), 
      "Coercing .+ to class [[:punct:]]character[[:punct:]]\\."
    )
    expect_equal(actual, expected)
  }
)

catz <- c(
  "catz/lolcat.gif",
  "moar cats/nyan cat.jpeg",
  "catz\\catz in loft\\ceiling cat.jpg",
  "catz/musical catz\\keyboard cat.bmp",
  "catbread.png",
  "kitties\\bonsai kitten.tiff",
  "kitties\\hipster kitty.pdf"
)

pwd <- getwd()
expected_catz <- create_expected_decomposed_path(
  dirname   = c(
    file.path(pwd, "catz"),
    file.path(pwd, "moar cats"),
    file.path(pwd, "catz/catz in loft"),
    file.path(pwd, "catz/musical catz"), getwd(),
    file.path(pwd, "kitties"),
    file.path(pwd, "kitties")
  ),
  filename  = c(
    "lolcat", "nyan cat", "ceiling cat", "keyboard cat", 
    "catbread", "bonsai kitten", "hipster kitty"
  ),
  extension = c(
    "gif", "jpeg", "jpg", "bmp", "png", "tiff", "pdf"
  ),
  row.names = catz
)

test_that(
  "decompose_path works with a character vector input.",
  {
    x <- catz
    expect_equal(decompose_path(x), expected_catz)
  }
)

test_that(
  "decompose_path works with a factor input.",
  {
    x <- factor(catz)
    expect_warning(
      actual <- decompose_path(x), 
      "Coercing .+ to class [[:punct:]]character[[:punct:]]\\."
    )
    expect_equal(actual, expected_catz)
  }
)
