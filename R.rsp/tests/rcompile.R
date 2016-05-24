library("R.rsp")
library("R.utils") # Arguments

verbose <- Arguments$getVerbose(TRUE)

path <- system.file(package="R.rsp")
path <- file.path(path, "rsp_tests")

pathname <- file.path(path, "trimming-1.txt.rsp")

verbose && enter(verbose, "rcompile()")

untils <- rev(eval(formals(R.rsp:::parse.RspParser)$until))
untils <- setdiff(untils, "*")

for (kk in seq_along(untils)) {
  until <- untils[kk]
  verbose && enter(verbose, sprintf("Until #%d ('%s') of %d", kk, until, length(untils)))
  s <- rcompile(file=pathname, until=until, output=RspString())
  verbose && ruler(verbose)
  cat(s)
  verbose && ruler(verbose)
  verbose && exit(verbose)
}

verbose && exit(verbose)
