library("R.rsp")
parse <- R.utils::parse # S3 generic parse() is exported by R.rsp
Arguments <- R.utils::Arguments
enter <- R.utils::enter
exit <- R.utils::exit

verbose <- Arguments$getVerbose(TRUE)

path <- system.file(package="R.rsp")
path <- file.path(path, "rsp_tests")

pathname <- file.path(path, "trimming-1.txt.rsp")

verbose && enter(verbose, "Validating that the RSP parse output can be deparsed")

untils <- rev(eval(formals(R.rsp:::parse.RspParser)$until))
untils <- setdiff(untils, "*")

for (kk in seq_along(untils)) {
  until <- untils[kk]
  verbose && enter(verbose, sprintf("Until #%d ('%s') of %d", kk, until, length(untils)))
  s0 <- rcompile(file=pathname, until=until, output=RspString())
  d0 <- parse(s0, until=until)
  stopifnot(inherits(d0, "RspDocument"))
  s1 <- asRspString(d0)
  d1 <- parse(s1, until=until)
  stopifnot(inherits(d1, "RspDocument"))
  stopifnot(identical(d1, d0))
  stopifnot(s1 == s0)
  verbose && exit(verbose)
}

verbose && exit(verbose)
