context("Registration repositories")

test_that("make_reg_url behaves", {
  expect_equal(make_reg_url("jefferislab/BridgingRegistrations"),
               "https://github.com/jefferislab/BridgingRegistrations")

  abc=c("a/b", "a/c")
  expect_equal(make_reg_url(abc), paste0("https://github.com/", abc))
})

test_that("add works",{
  dir.create(regroot<-tempfile())
  on.exit(unlink(regroot, recursive = TRUE))
  op=options(nat.templatebrains.regdirs='testing')
  on.exit(options(op), add = TRUE)
  add_reg_folders(regroot)
  rop=getOption('nat.templatebrains.regdirs')
  expect_equal(length(rop), 2L)
  expect_equal(rop[2], "testing")
  add_reg_folders(regroot)
  expect_equal(getOption('nat.templatebrains.regdirs'), rop)
  options(nat.templatebrains.regdirs='testing')
  add_reg_folders(regroot, first = FALSE)
  expect_equal(getOption('nat.templatebrains.regdirs'), rev(rop))

  expect_null(add_reg_folders(character(0)))
})

test_that("cloning registrations works", {
  skip_on_cran()
  skip_if_not_installed('git2r')
  op=options(nat.templatebrains.regdirs=NULL)
  on.exit(options(op), add = TRUE)
  # remove it in case it was already there!
  unlink(local_reg_dir_for_url('https://github.com/jefferislab/TestRegRepo'), recursive = TRUE)
  download_reg_repo("jefferislab/TestRegRepo")
  expect_equal(getOption('nat.templatebrains.regdirs'),
               local_reg_dir_for_url('https://github.com/jefferislab/TestRegRepo'))
  expect_true(update_result<-update_reg_repos()@up_to_date)
  unlink(local_reg_dir_for_url('https://github.com/jefferislab/TestRegRepo'), recursive = TRUE)
})
