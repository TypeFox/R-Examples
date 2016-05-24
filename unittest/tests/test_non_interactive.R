R_binary <- file.path(R.home("bin"), "R")

run_script <- function(script, expected_status, expected_out, description) {
    # solaris does not like pipes so use tmp files as intermediaries
    tmpfiles <- tempfile(pattern = c('R_unittest_stdout_','R_unittest_stderr_'), tmpdir = tempdir())
    exit_status <- withCallingHandlers(
        system2(
            R_binary,
            c("--vanilla", "--slave"),
            input=script,
            wait = TRUE, stdout = tmpfiles[1], stderr = tmpfiles[2]),
        warning = function (w) {
            invokeRestart("muffleWarning")
        }
    )
    actual <- readLines(tmpfiles[1])  # only interested in stdout
    if( isTRUE(all.equal(actual, expected_out)) && exit_status == expected_status) {
        cat("ok\n")
    } else {
        cat("\nExpected status",
            expected_status,
            "\nGot status",
            exit_status,
            "\nExpected stdout:",
            expected_out,
            "\nGot stdout:",
            actual,
            sep = "\n"
        )
        stop( description )
    }
    invisible(c(exit_status, actual))
}

# one test one success
run_script(
    "library(unittest, quietly = TRUE)\nok(1==1,\"1 equals 1\")",
    0,
    c(
        "ok - 1 equals 1",
        "# Looks like you passed all 1 tests."
    ),
    "One test one success case not as expected"
) 

# Success with a multi-line expression
run_script(
    "library(unittest, quietly = TRUE)\nok(all.equal(c('This is a string', 'This is a string too', 'Exciting times'),\nc('This is a string', 'This is a string too', 'Exciting times')))",
    0,
    c(
        'ok - all.equal(c("This is a string", "This is a string too", "Exc',
        "# Looks like you passed all 1 tests."
    ),
    "One test one success case not as expected"
) 

# two tests two sucesses
run_script(
    "library(unittest, quietly = TRUE)\nok(1==1,\"1 equals 1\")\nok(2==2,\"2 equals 2\")",
    0,
    c(
        "ok - 1 equals 1",
        "ok - 2 equals 2",
        "# Looks like you passed all 2 tests."
    ),
    "Two tests two successes case not as expected"
) 

# one test one failure 
run_script(
    "library(unittest, quietly = TRUE)\nok(1!=1,\"1 equals 1\")",
    10,
    c(
        "not ok - 1 equals 1",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "# Looks like you failed 1 of 1 tests."
    ),
    "One test one failure case not as expected"
)

# four tests two failures 
run_script(
    "library(unittest, quietly = TRUE)\nok(1==1,\"1 equals 1\")\nok(2!=2,\"2 equals 2\")\nok(3==3,\"3 equals 3\")\nok(4!=4,\"4 equals 4\")",
    10,
    c(
        "ok - 1 equals 1",
        "not ok - 2 equals 2",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "ok - 3 equals 3",
        "not ok - 4 equals 4",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "# Looks like you failed 2 of 4 tests."
    ),
    "Four tests two failures case not as expected"
)

# check detaching stops non_interactive_exit functionality
run_script(
    "library(unittest, quietly = TRUE)\nok(1!=1,\"1 equals 1\")\ndetach(package:unittest,unload=FALSE)",
    0,
    c(
        "not ok - 1 equals 1",
        "# Test returned non-TRUE value:",
        "# [1] FALSE"
    ),
    "detaching stops non_interactive_exit functionality"
)

# and if we re-attach it works again
run_script(
    "library(unittest, quietly = TRUE)\nok(1!=1,\"1 equals 1\")\ndetach(package:unittest,unload=FALSE)\nlibrary(unittest, quietly = TRUE)\nok(2!=2,\"2 equals 2\")",
    10,
    c(
        "not ok - 1 equals 1",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "not ok - 2 equals 2",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "# Looks like you failed 1 of 1 tests."
    ),
    "detaching stops non_interactive_exit functionality and then re-attaching resets and the rest still works"
)

# check detaching and unloading stops non_interactive_exit functionality
run_script(
    "library(unittest, quietly = TRUE)\nok(1!=1,\"1 equals 1\")\ndetach(package:unittest,unload=TRUE)",
    0,
    c(
        "not ok - 1 equals 1",
        "# Test returned non-TRUE value:",
        "# [1] FALSE"
    ),
    "detaching and unloading stops non_interactive_exit functionality"
)

# and if we reload and re-attach it works again
run_script(
    "library(unittest, quietly = TRUE)\nok(1!=1,\"1 equals 1\")\ndetach(package:unittest,unload=TRUE)\nlibrary(unittest, quietly = TRUE)\nok(2!=2,\"2 equals 2\")",
    10,
    c(
        "not ok - 1 equals 1",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "not ok - 2 equals 2",
        "# Test returned non-TRUE value:",
        "# [1] FALSE",
        "# Looks like you failed 1 of 1 tests."
    ),
    "detaching and unloading stops non_interactive_exit functionality and then reloading and re-attaching resets and the rest still works"
)

# Failure outside test
# NB: The error message is on stderr, not stdout, so the TAP output lies. Is this bad?
run_script(
    paste(
        "library(unittest, quietly = TRUE)",
        "ok(1==1, '1 equals 1')",
        "stop('eek\nook')",
        "ok(2==2, '2 equals 2')",
        "", sep = "\n"
    ),
    1,
    c(
        "ok - 1 equals 1",
        "# Looks like you passed all 1 tests.",
        NULL
    ),
    "Failure outside tests"
)

# tryCatch() doesn't count as failure
run_script(
    paste(
        "library(unittest, quietly = TRUE)",
        "ok(1==1, '1 equals 1')",
        "tryCatch(stop('not fatal'), error = function (e) NULL)",
        "ok(2==2, '2 equals 2')",
        "", sep = "\n"
    ),
    0,
    c(
        "ok - 1 equals 1",
        "NULL",
        "ok - 2 equals 2",
        "# Looks like you passed all 2 tests.",
        NULL
    ),
    "Caught errors outside tests"
)
