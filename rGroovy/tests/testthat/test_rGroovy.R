#
# IF Error: '\.' is an unrecognized escape in character string starting "'\."
#    AND the problem involves build -> test package
#
# THEN TRY THE FOLLOWING:
#
#     install.packages("devtools")
#
# TODO:
#   - install.packages("testthat")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# WARNING:
#
# Using Test Package will result in a java.lang.ClassNotFoundException -- to run this test use Check Package.
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#

#
# We do not want these tests being run on CRAN because the Groovy dependencies must be included otherwise they will
# fail. We can test for this by uncommenting the following line and when devtools::check is executed, we should not see
# any failed tests.
#
# NOTE: This code should only be used when testing locally.
#
# Sys.setenv(NOT_CRAN='false')

test_that (
    "Calling Initialize twice does not raise an error.",
    {
        skip_on_cran()

        Initialize ()
        expect_warning (Initialize ())
    }
)

test_that (
    "The Groovy script executes correctly.",
    {
        skip_on_cran()

        groovyScript <- "return 'Hello world!'"

        result <- Evaluate (groovyScript=groovyScript)

        expect_equal("Hello world!", result)
    }
)

test_that (
    "A NULL Groovy script raises an exception.",
    {
        skip_on_cran()

        expect_error(Evaluate (groovyScript=NULL), "The groovyScript parameter cannot be NULL.")
    }
)

test_that (
    "Execute with NULL variables does not raise an exception.",
    {
        skip_on_cran()

        groovyScript <- "return 'Hello world!'"

        Execute (groovyScript=groovyScript, variables=NULL)
    }
)