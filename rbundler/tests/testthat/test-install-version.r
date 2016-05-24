context('rbundler can install a specific version of a package')

repos <- 'http://cran.rstudio.com'
options(repos=repos)
options(pkgType = 'source')

dependency <- mock_dependency(repos=repos)

create_new_lib <- function() {
  lib <- file.path(tempdir(), 'library', as.numeric(Sys.time()))

  dir.create(lib, recursive=TRUE)
  .libPaths(lib)
  lib
}

test_that('can compare an available version with a requested version', {
  expect_true(compare_versions('1', '==', '1'))
  expect_false(compare_versions('2', '==', '1'))
  expect_true(compare_versions('2', '>=', '1'))
  expect_true(compare_versions('0.2.2', '>=', '0.2.1'))
})

test_that('the correct version for installation can be determined',{

  old_versions <- c('1', '2')
  available_versions <- c(old_versions, '3')
  expect_version_to_install <- function(expected, version, comparator) {
    version_to_install <- determine_version_to_install(available_versions, version, comparator)
    expect_equal(
      version_to_install,
      expected,
      info = sprintf('Expected version [%s] but got [%s]', expected, paste(version_to_install, collapse=','))
    )
  }


  for(version in available_versions) {
    expect_version_to_install('3', version, NA)
  }

  for(version in available_versions) {
    expect_version_to_install(version, version, '==')
    expect_version_to_install(version, version, '<=')
  }

  for(version in old_versions) {
    expect_version_to_install('3', version, '>=')
    expect_version_to_install('3', version, '>')
  }

  expect_version_to_install('1', '2', '<')
  expect_version_to_install('2', '3', '<')

})

test_that('can validate compare clause', {
  expect_that({validate_compare(NULL)}, throws_error())
})

test_that('can read archive from both cran and non-cran-like repsitories', {
  read_archive_rds('http://cran.rstudio.com') # Check that there are no errors
  expect_equal(read_archive_rds('http://cran.does.not.exist'), list())
})

test_that('can install an old package version', {

  lib <- create_new_lib()
  install_version(dependency$name, dependency$version, '==')

  expected_path <- file.path(lib, dependency$name)
  package <- as.package(expected_path)

  expect_true(file.exists(expected_path))
  expect_equal(package$version, dependency$version)

})

test_that('install_version is idempotent', {

  create_new_lib()
  # Try to install the package several times, each time with a compatible version:
  expect_true(install_version(dependency$name, dependency$version, '=='))
  expect_false(install_version(dependency$name, dependency$version, '=='))
  expect_false(install_version(dependency$name, dependency$version, '>='))
  expect_false(install_version(dependency$name, dependency$version, '<='))

})

test_that('install_version throws an exception if a package has already been installed but has the wrong version', {

  create_new_lib()
  # Try to install the package several times. The first time with a valid version, and subsequent times with invalid
  # compatibility:
  expect_true(install_version(dependency$name, dependency$version, '=='))
  expect_error(install_version(dependency$name, dependency$version, '<'))
  expect_error(install_version(dependency$name, dependency$version, '>'))

})

test_that('correct package is returned for mac, windows, and source package versions', {

  library('tools')

  expect_package_extension <- function(type, expected_extension) {
    available_versions <- find_available_versions(dependency$name, type=type)
    extension <- file_ext(available_versions[available_versions$source == 'root', ]$url)
    expect_equal(extension, expected_extension)
  }

  expect_package_extension(type='source', expected_extension='gz')
  expect_package_extension(type='mac.binary', expected_extension='tgz')
  expect_package_extension(type='win.binary', expected_extension='zip')

}) 

