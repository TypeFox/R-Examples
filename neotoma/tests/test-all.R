#  Code as used in rOpenSci package rfigshare
#  author: Carl Boettiger
#  file: https://github.com/ropensci/rfigshare/commit/4f479b6e34f2ab5e89b3aa3da2355156ae659950#diff-aec6d055aaf98780cdcb2bb54dd3631c
#  minor modification by : S. Goring

if (packageVersion("testthat") >= "0.7.1.99") {
  library("testthat")
  test_check("neotoma")
} 
if (packageVersion("testthat") < "0.7.1.99") {
  library("testthat")
  test_package("neotoma")
}
