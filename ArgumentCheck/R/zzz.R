.onAttach <- function(libname, pkgname)
{
  packageStartupMessage("ArgumentCheck has been deprecated and is no longer ",
                        "actively maintained.\n",
                        "I recommend switching to the `checkmate` package:\n",
                        "https://cran.r-project.org/web/packages/checkmate/index.html\n\n",
                        "Should you need help migrating your code, please contact me at \n",
                        "benjamin.nutter@gmail.com\n",
                        "\n",
                        "ArgumentCheck will be removed from CRAN in early 2017.")
}