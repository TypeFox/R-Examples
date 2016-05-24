PACKAGE_NAME <- methods::getPackageName()
PACKAGE_VERSION <- "0.1"

auto_extdata <-
  function (assign.env = parent.frame())
  {
    extension_pattern <- "[.]rds$"
    files <- dir(extdata_path(assign.env), pattern = extension_pattern,
                 full.names = TRUE)
    read_rds(files, assign.env = assign.env)
  }

read_rds <-
  function (..., assign.env = parent.frame())
  {
    dots <- c(...)
    if (length(dots) == 0L) {
      return()
    }
    unnamed <- (names(dots) == "")
    if (length(unnamed) == 0L) {
      unnamed <- rep(TRUE, length(dots))
    }
    names(dots)[unnamed] <- name_from_rds(dots[unnamed])
    dots_expr <- lapply(dots, function(dot) {
      lazyeval::lazy_(bquote({
        f <- tempfile("rextdata", fileext = ".rds")
        on.exit(unlink(f))
        curl::curl_download(.(dot), f)
        readRDS(f)
      }), baseenv())
    })
    delayed_assign_(.dots = dots_expr, assign.env = assign.env)
  }

extdata_path <-
  function (package.env)
  {
    system.file(extdata_name(), package = utils::packageName(package.env))[[1L]]
  }

extdata_name <-
  function ()
    c("extdata", file.path("inst", "extdata"))

name_from_rds <-
  function (x)
  {
    if (length(x) == 0L)
      return()
    gsub("[.]rds$", "", basename(x))
  }

delayed_assign_ <-
  function (..., .dots, assign.env = parent.frame())
  {
    dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
    ret <- mapply(delayed_assign_one, names(dots), dots, MoreArgs = list(assign.env = assign.env))
    invisible(ret)
  }

delayed_assign_one <-
  function (name, expr, assign.env = parent.frame())
  {
    force(name)
    force(expr)
    delayedAssign(name, lazyeval::lazy_eval(expr), assign.env = assign.env)
    invisible(name)
  }

read_rds(
  file.path(
    paste0("https://raw.githubusercontent.com/krlmlr/", PACKAGE_NAME, "/v",
           PACKAGE_VERSION, "/inst/extdata"),
    c("break_even.rds", "p_values_7.rds", "p_values_agg_agg.rds",
      "p_values_agg.rds", "timings.rds")
  )
)
