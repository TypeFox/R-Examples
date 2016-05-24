# Unit test: "stratify"
# Tests error handling
test_stratify <-
function() {

    # load rgdal
    if (suppressWarnings(!require(rgdal))) {
        stop("This unit test requires package 'rgdal'.\nThis package is currently not available. Please install 'rgdal' first.", call. = FALSE)
    }    
    
    # read vector map for testing
    map <- readOGR(dsn = system.file("maps", package = "spcosa"), layer = "farmsum")

    # define prior points
    priorPoints <- data.frame(x = 1:10, y = 1:10)
    coordinates(priorPoints) <- ~ x * y

    # check exception: nStrata > nGridCells
    checkException(
        stratify(
            object = map,
            nStrata = 40,
            priorPoints = NULL,
            maxIterations = 1000,
            nTry = 1,
            nGridCells = 5,
            equalArea = FALSE
        )
    )

    # check exception: nStrata is not a scalar
    checkException(
        stratify(
            object = map,
            nStrata = 1:5,
            priorPoints = NULL,
            maxIterations = 1000,
            nTry = 1,
            nGridCells = 10000,
            equalArea = FALSE
        )
    )

    # check exception: nStrata < 1
    checkException(
        stratify(
            object = map,
            nStrata = 0,
            priorPoints = NULL,
            maxIterations = 1000,
            nTry = 1,
            nGridCells = 10000,
            equalArea = FALSE
        )
    )

    # check exception: priorPoints is not an instance of "SpatialPoints"
    checkException(
        stratify(
            object = map,
            nStrata = 10,
            priorPoints = matrix(runif(10), nrow = 5, ncol = 2),
            maxIterations = 1000,
            nTry = 1,
            nGridCells = 10000,
            equalArea = FALSE
        )
    )

    # check exception: test should fail for priorPoints in combination of equalArea = TRUE
    checkException(
        stratify(
            object = map,
            nStrata = 10,
            priorPoints = priorPoints,
            maxIterations = 1000,
            nTry = 1,
            nGridCells = 10000,
            equalArea = TRUE
        )
    )

    # check exception: test should fail in the rare case that 'priorPoints' does not have coordinates
    checkException(
        stratify(
            object = map,
            nStrata = 10,
            priorPoints = suppressWarnings(priorPoints[-seq_len(nrow(coordinates(priorPoints))), ]),
            maxIterations = 1000,
            nTry = 1,
            nGridCells = 10000,
            equalArea = FALSE
        )
    )
}
