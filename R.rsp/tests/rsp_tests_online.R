# Find all RSP files
path <- system.file("rsp_tests_online", package="R.rsp")
pathnames <- list.files(path=path, pattern="[.]R$", full.names=TRUE)
pathnames <- pathnames[file_test("-f", pathnames)]

for (kk in seq_along(pathnames)) {
  source(pathnames[kk], echo=TRUE)
}
