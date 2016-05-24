context("A developer can bundle a package and it's dependencies.")

repos <- 'http://cran.rstudio.com'
options(repos = repos)
options(pkgType = 'source')

#' Ensures that a package and it's dependencies are installed when running 'bundle'
#' @param desc a description of the test
#' @param package the package to test
#' @param expected_dependencies A vector of dependencies which must be installed
test_bundle <- function(desc, package, expected_dependencies, confirm_dependency_versions = FALSE) {

  test_that(
            desc,
            {

              lib <- file.path(package$path, ".Rbundle")

              bundle_and_validate <- function(overwrite) {

                bundle(pkg = package$path, bundle_path = lib, overwrite = overwrite)

                expect_true(
                            file.exists(lib),
                            info=sprintf("Bundler library [%s] was not created.", lib)
                            )

                bundle_package_path <- file.path(lib, package$package)

                expect_true(
                            file.exists(bundle_package_path),
                            sprintf(
                                    "Package was not successfully installed into bundler library: [%s]",
                                    bundle_package_path
                                    )
                            )

                for(dependency in expected_dependencies) {

                  bundle_dependency_path <- file.path(lib, dependency$name)

                  expect_true(
                              file.exists(bundle_dependency_path),
                              sprintf(
                                      "Dependency was not successfully installed into bundler library: [%s]",
                                      bundle_dependency_path
                                      )
                              )

                  if(confirm_dependency_versions) {

                    package <- as.package(bundle_dependency_path)
                    expect_equal(package$version, dependency$version)

                  }

                }

                # Note: we're checking that only the basename of the library is in libPaths, due to an issue with OS X
                # having '/var' point to '/private/var'. R thinks it's using 'var', when in reality, .libPaths
                # points to '/private/var'.
                expect_true(basename(lib) %in% basename(.libPaths()), info = sprintf("Did not find [%s] in .libPaths().", lib))

              }

              # First bundle should install all packages and create the .Rbundle directory.
              bundle_and_validate(overwrite = FALSE)
              # Second bundle - with no overwrite - should install nothing, and effectively be idempotent
              bundle_and_validate(overwrite = FALSE)
              # Third bundle - this time with overwrite - should delete the bundle library and re-install everything
              bundle_and_validate(overwrite = TRUE)

              # Reset libraries for the next test.
              .libPaths('new')

            }
            )
}

test_path <- file.path(tempdir(), sprintf('bundle-test-%s', as.numeric(Sys.time())))
dir.create(test_path)
dependency <- mock_dependency(repos=repos)
mock_packages <- create_mock_packages(test_path, dependency, repos)

test_bundle(
            desc = "Bundling a package with no dependencies creates a bundle directory with the project installed and no dependencies installed.",
            package = mock_packages[['nodependencies']],
            expected_dependencies = list()
            )

test_bundle(
            desc = "Bundling a package with dependencies creates a bundle directory with the project installed and all dependencies installed.",
            package = mock_packages[["simpledependencies"]],
            expected_dependencies=list(dependency)
            )

test_bundle(
            desc = "Bundling a package with versioned dependencies creates a bundle directory with the project installed, all dependencies installed, and correct dependency versions.",
            package = mock_packages[["versioneddependencies"]],
            expected_dependencies=list(dependency),
            confirm_dependency_versions = TRUE
            )

