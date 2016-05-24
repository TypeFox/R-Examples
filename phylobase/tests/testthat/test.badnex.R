#
# --- Test badnex.R ---
#

test_that("Malformed Nexus File should not work.", {    
    if (Sys.getenv("RCMDCHECK") == FALSE) {
        pth <- file.path(getwd(), "..", "inst", "nexusfiles")
    } else {
        pth <- system.file(package="phylobase", "nexusfiles")
    }
    badFile <- file.path(pth, "badnex.nex")
    expect_error(readNexus(file=badFile))
})


