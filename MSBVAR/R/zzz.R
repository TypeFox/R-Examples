# Start up functions for the MSBVAR package
#
# 2012-01-19 : updated to use new loading mechanism.
# 2012-05-20 : automated build date
# 2014-03-16 : New updates to LICENSE and version
# 2014-06-09 : Updated to use packageStartupMessage(); allows user
#              suppression if required


.onAttach <- function(...)
{
    dt <- date()
    x <- regexpr("[0-9]{4}", dt)
    yr <- substr(dt, x[1], x[1] + attr(x, "match.length") - 1)

    packageStartupMessage(
        c("##\n## MSBVAR Package v.0.9-2\n",
          paste("## Build date: ", date(), "\n"),
          paste("## Copyright (C) 2005-", yr, ", Patrick T. Brandt\n", sep=""),
          "## Written by Patrick T. Brandt\n",
          "##\n## Support provided by the U.S. National Science Foundation\n",
          "## (Grants SES-0351179, SES-0351205, SES-0540816, and SES-0921051)\n##\n"
          ))

    rm(dt, x, yr)
}

.onUnload <- function(libpath) {
    library.dynam.unload("MSBVAR", libpath)
}
