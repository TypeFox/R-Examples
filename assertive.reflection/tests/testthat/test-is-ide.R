# TODO: tests for is_rstudio_current/desktop/server

test_that("test.is_architect.some_ide.returns_true_if_ide_is_architect",
{
  rj_is_loaded <- "package:rj" %in% search()
  device_name <- formals(getOption("device"))$name
  expected <- rj_is_loaded && !is.null(device_name) && device_name == "rj.gd"
  actual <- is_architect()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), noquote("You are not running Architect/StatET."))
  }  
})

test_that("test.is_revo_r.any_os.returns_true_if_ide_is_revo_r", {
  expected <- exists("Revo.version", "package:base", inherits = FALSE) &&
    is.list(get("Revo.version", "package:base", inherits = FALSE))
  actual <- is_revo_r()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), noquote("You are not running Revolution R."))
  }
})

test_that("test.is_rstudio.any_os.returns_true_if_ide_is_rstudio", 
{
  env <- Sys.getenv("RSTUDIO", "0")
  expected <- env == "1"
  actual <- is_rstudio()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(cause(actual), noquote("You are not running RStudio."))
  }
})

test_that(
  "test.is_rstudio_desktop.any_os.returns_true_if_ide_is_rstudio", 
  {
    gui <- .Platform$GUI
    actual <- is_rstudio_desktop()
    if(!is_rstudio()) 
    {
      expect_false(strip_attributes(actual))
      expect_equal(cause(actual), noquote("You are not running RStudio."))
    } else if("tools:rstudio" %in% search())
    {
      e <- as.environment("tools:rstudio")
      if(!".rs.api.versionInfo" %in% ls(e, all.names = TRUE))
      {
        expect_equal(strip_attributes(actual), NA)
        expect_equal(
          cause(actual), 
          noquote("You are using an old version of RStudio, which does not tell you version information.")
        )
      } else
      {
        info <- e$.rs.api.versionInfo()
        expect_equal(strip_attributes(actual), info$mode == "desktop")
        if(!actual)
        {
          expect_equal(
            cause(actual), 
            noquote("You are running the server version of RStudio.")
          )
        }
      }
    } else
    {
      expect_equal(strip_attributes(actual), NA)
      expect_equal(
        cause(actual), 
        noquote("'tools:rstudio' is not loaded, so the RStudio API is not available.")
      )
    }
  }
)

test_that(
  "test.is_rstudio_server.any_os.returns_true_if_ide_is_rstudio", 
  {
    gui <- .Platform$GUI
    actual <- is_rstudio_server()
    if(!is_rstudio()) 
    {
      expect_false(strip_attributes(actual))
      expect_equal(cause(actual), noquote("You are not running RStudio."))
    } else if("tools:rstudio" %in% search())
    {
      e <- as.environment("tools:rstudio")
      if(!".rs.api.versionInfo" %in% ls(e, all.names = TRUE))
      {
        expect_equal(strip_attributes(actual), NA)
        expect_equal(
          cause(actual), 
          noquote("You are using an old version of RStudio, which does not tell you version information.")
        )
      } else
      {
        info <- e$.rs.api.versionInfo()
        expect_equal(strip_attributes(actual), info$mode == "server")
        if(!actual)
        {
          expect_equal(
            cause(actual), 
            noquote("You are running the desktop version of RStudio.")
          )
        }
      }
    } else
    {
      expect_equal(strip_attributes(actual), NA)
      expect_equal(
        cause(actual), 
        noquote("'tools:rstudio' is not loaded, so the RStudio API is not available.")
      )
    }
  }
)




