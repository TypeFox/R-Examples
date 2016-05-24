context("Pooling")



test_that("no pooling occurs if all sample areas are within the target area.", {

        targetArea <- c(0.09, 0.11)
        nSamples <- 10L
        area <- rep.int(x = 0.10, times = nSamples)
        sampleId <- 1:nSamples
    
    	expect_that(
            pool(sampleId = sampleId, area = area, targetArea = targetArea),
            shows_message("No pooling necessary.")
        )
        expect_that(
            pool(sampleId = sampleId, area = area, targetArea = targetArea),
            is_identical_to(sampleId)
        )
})



test_that("no pooling occurs if all sample areas are greater than the target area.", {

        targetArea <- c(0.09, 0.11)
        nSamples <- 10L
        area <- rep.int(x = 0.20, times = nSamples)
        sampleId <- 1:nSamples
    
        expect_that(
            pool(sampleId = sampleId, area = area, targetArea = targetArea),
            throws_error("No pooling possible")
        )
})



test_that("pooling occurs for trivial cases.", {

        targetArea <- c(0.09, 0.11)
        nSamples <- 10L
        area <- rep.int(x = 0.01, times = nSamples)
        sampleId <- 1:nSamples
    
        expect_that(
            pool(sampleId = sampleId, area = area, targetArea = targetArea),
            is_identical_to(rep.int(x = 1L, times = nSamples))
        )
})



test_that("no pooling occurs if the sum of all sample areas is below the target area.", {

        targetArea <- c(0.09, 0.11)
        nSamples <- 10L
        area <- rep.int(x = 0.001, times = nSamples)
        sampleId <- 1:nSamples
    
        expect_that(
            pool(sampleId = sampleId, area = area, targetArea = targetArea),
            throws_error("No pooling possible")
        )
})



test_that("no pooling occurs if all sample areas are slightly smaller than the target area.", {

        targetArea <- c(0.09, 0.11)
        nSamples <- 10L
        area <- rep.int(x = 0.08, times = nSamples)
        sampleId <- 1:nSamples
    
        expect_that(
            pool(sampleId = sampleId, area = area, targetArea = targetArea),
            is_identical_to(rep.int(x = NA_integer_, times = nSamples))
        )
})



test_that("pooling processes all samples if possible.", {

        targetArea <- c(0.09, 0.11)
        nSamples <- 8L
        area <- rep.int(x = 0.025, times = nSamples)
        sampleId <- 1:nSamples
    
        for (i in 1:10) {
            expect_that(
                any(is.na(pool(sampleId = sampleId, area = area, targetArea = targetArea))),
                is_false()
            )
        }
})



test_that("areas of pools are in target interval.", {
        targetArea <- c(0.09, 0.11)
        for (i in 1:10) {
            nSamples <- sample(x = 10:250, size = 1)
            sampleId <- 1:nSamples
            area <- runif(n = nSamples, min = 0.01, max = 0.04)
            index <- pool(sampleId = sampleId, area = area, targetArea = targetArea)
            expect_that(
                all(tapply(X = area, INDEX = index, FUN = sum) %inInterval% targetArea),
                is_true()
            )
        }
})