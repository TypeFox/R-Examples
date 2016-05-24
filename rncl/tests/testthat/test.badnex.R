#
# --- Test badnex.R ---
#

test_that("Malformed Nexus File should not work.", {
    if (Sys.getenv("RCMDCHECK") == FALSE) {
        pth <- file.path(getwd(), "..", "inst", "nexusfiles")
    } else {
        pth <- system.file(package="rncl", "nexusfiles")
    }
    badFile <- file.path(pth, "badnex.nex")
    stopifnot(file.exists(badFile))
    expect_error(make_phylo(file=badFile))
})
