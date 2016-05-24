context("BEDMatrix")

parseRaw <- function(path) {
    file <- file(path)
    lines <- readLines(file)
    out <- t(sapply(strsplit(lines[2:length(lines)], " "), function(line) {
        line <- line[7:length(line)]
        line <- gsub("NA", NA, line)
        as.integer(line)
    }))
    rownames(out) <- paste0("id_", 1:nrow(out))
    colnames(out) <- paste0("mrk_", 1:ncol(out))
    close(file)
    return(out)
}

# Prepare dummy data
raw <- parseRaw(system.file("extdata", "example.raw", package = "BEDMatrix"))
examplePath <- system.file("extdata", "example.bed", package = "BEDMatrix")
standalonePath <- system.file("extdata", "standalone.bed", package = "BEDMatrix")

for (path in c(examplePath, sub(".bed", "", examplePath))) {

    test_that("it throws an error if file does not exist", {
        expect_error(BEDMatrix("NOT_FOUND"), "File not found\\.")
    })

    test_that("it throws an error if file is not a BED file", {
        expect_error(BEDMatrix(system.file("extdata", "example.raw", package = "BEDMatrix")), "File is not a binary PED file\\.")
    })

    test_that("it determines n from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(nrow(bed), 6)
        expect_message(BEDMatrix(path = path), "Extracting number of individuals and rownames from FAM file\\.\\.\\.")
    })

    test_that("it throws an error if FAM file is not found and n is not given", {
        expect_error(BEDMatrix(path = standalonePath), "FAM file of same name not found\\. Provide number of individuals \\(n\\)\\.")
    })

    test_that("it determines rownames from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(rownames(bed), c("1_1", "1_2", "1_3", "2_1", "2_2", "2_3"))
        expect_message(BEDMatrix(path = path), "Extracting number of individuals and rownames from FAM file\\.\\.\\.")
    })

    test_that("it determines p from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(ncol(bed), 3)
        expect_message(BEDMatrix(path = path), "Extracting number of markers and colnames from BIM file\\.\\.\\.")
    })

    test_that("it throws an error if BIM file is not found and p is not given", {
        expect_error(BEDMatrix(path = standalonePath), "FAM file of same name not found\\. Provide number of individuals \\(n\\)\\.")
    })

    test_that("it determines colnames from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(colnames(bed), c("snp1_G", "snp2_1", "snp3_A"))
        expect_message(BEDMatrix(path = path), "Extracting number of markers and colnames from BIM file\\.\\.\\.")
    })

    test_that("it accepts n and p if FAM or BIM file are present", {
        bed <- BEDMatrix(path = path, n = 6, p = 3)
        expect_equal(dimnames(bed), list(NULL, NULL))
    })

    test_that("it accepts n and p if FAM or BIM file is not found", {
        bed <- BEDMatrix(path = standalonePath, n = 6, p = 3)
        expect_equal(dimnames(bed), list(NULL, NULL))
    })

    test_that("it throws an error if dimensions are wrong", {
        expect_error(BEDMatrix(path = path, n = 10, p = 5), "n or p does not match the dimensions of the file\\.")
    })

}

# Prepare dummy BED matrix
bed <- BEDMatrix(path = examplePath, n = 6, p = 3)
rownames(bed) <- paste0("id_", 1:nrow(bed))
colnames(bed) <- paste0("mrk_", 1:ncol(bed))

test_that("subsetting", {

    expect_equal(bed[], raw[])

    expect_equal(bed[1, ], raw[1, ])
    expect_equal(bed[, 1], raw[, 1])
    expect_equal(bed[1, 1], raw[1, 1])
    expect_equal(bed[1, , drop = FALSE], raw[1, , drop = FALSE])
    expect_equal(bed[, 1, drop = FALSE], raw[, 1, drop = FALSE])
    expect_equal(bed[1, 1, drop = FALSE], raw[1, 1, drop = FALSE])

    expect_equal(bed[1:2, ], raw[1:2, ])
    expect_equal(bed[, 1:2], raw[, 1:2])
    expect_equal(bed[1:2, 1:2], raw[1:2, 1:2])
    expect_equal(bed[1:2, , drop = FALSE], raw[1:2, , drop = FALSE])
    expect_equal(bed[, 1:2, drop = FALSE], raw[, 1:2, drop = FALSE])
    expect_equal(bed[1:2, 1:2, drop = FALSE], raw[1:2, 1:2, drop = FALSE])

    expect_equal(bed[2:1, ], raw[2:1, ])
    expect_equal(bed[, 2:1], raw[, 2:1])
    expect_equal(bed[2:1, 2:1], raw[2:1, 2:1])
    expect_equal(bed[2:1, , drop = FALSE], raw[2:1, , drop = FALSE])
    expect_equal(bed[, 2:1, drop = FALSE], raw[, 2:1, drop = FALSE])
    expect_equal(bed[2:1, 2:1, drop = FALSE], raw[2:1, 2:1, drop = FALSE])

    expect_equal(bed[c(3, 1), ], raw[c(3, 1), ])
    expect_equal(bed[, c(3, 1)], raw[, c(3, 1)])
    expect_equal(bed[c(3, 1), c(3, 1)], raw[c(3, 1), c(3, 1)])
    expect_equal(bed[c(3, 1), , drop = FALSE], raw[c(3, 1), , drop = FALSE])
    expect_equal(bed[, c(3, 1), drop = FALSE], raw[, c(3, 1), drop = FALSE])
    expect_equal(bed[c(3, 1), c(3, 1), drop = FALSE], raw[c(3, 1), c(3, 1), drop = FALSE])

    expect_equal(bed[c(TRUE, FALSE, TRUE), ], raw[c(TRUE, FALSE, TRUE), ])
    expect_equal(bed[, c(TRUE, FALSE, TRUE)], raw[, c(TRUE, FALSE, TRUE)])
    expect_equal(bed[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], raw[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)])
    expect_equal(bed[c(TRUE, FALSE, TRUE), , drop = FALSE], raw[c(TRUE, FALSE, TRUE), , drop = FALSE])
    expect_equal(bed[, c(TRUE, FALSE, TRUE), drop = FALSE], raw[, c(TRUE, FALSE, TRUE), drop = FALSE])
    expect_equal(bed[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], raw[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE])

    expect_equal(bed[1], raw[1])
    expect_equal(bed[1:2], raw[1:2])
    expect_equal(bed[2:1], raw[2:1])
    expect_equal(bed[c(3, 1)], raw[c(3, 1)])
    expect_equal(bed[c(TRUE, FALSE, TRUE)], raw[c(TRUE, FALSE, TRUE)])
    expect_equal(bed[raw > 1], raw[raw > 1])

    expect_equal(bed["id_1", ], raw["id_1", ])
    expect_equal(bed[, "mrk_1"], raw[, "mrk_1"])
    expect_equal(bed["id_1", "mrk_1"], raw["id_1", "mrk_1"])
    expect_equal(bed["id_1", , drop = FALSE], raw["id_1", , drop = FALSE])
    expect_equal(bed[, "mrk_1", drop = FALSE], raw[, "mrk_1", drop = FALSE])
    expect_equal(bed["id_1", "mrk_1", drop = FALSE], raw["id_1", "mrk_1", drop = FALSE])

    expect_equal(bed[c("id_1", "id_2"), ], raw[c("id_1", "id_2"), ])
    expect_equal(bed[, c("mrk_1", "mrk_2")], raw[, c("mrk_1", "mrk_2")])
    expect_equal(bed[c("id_1", "id_2"), c("mrk_1", "mrk_2")], raw[c("id_1", "id_2"), c("mrk_1", "mrk_2")])
    expect_equal(bed[c("id_1", "id_2"), , drop = FALSE], raw[c("id_1", "id_2"), , drop = FALSE])
    expect_equal(bed[, c("mrk_1", "mrk_2"), drop = FALSE], raw[, c("mrk_1", "mrk_2"), drop = FALSE])
    expect_equal(bed[c("id_1", "id_2"), c("mrk_1", "mrk_2"), drop = FALSE], raw[c("id_1", "id_2"), c("mrk_1", "mrk_2"), drop = FALSE])

    expect_equal(bed[c("id_2", "id_1"), ], raw[c("id_2", "id_1"), ])
    expect_equal(bed[, c("mrk_2", "mrk_1")], raw[, c("mrk_2", "mrk_1")])
    expect_equal(bed[c("id_2", "id_1"), c("mrk_2", "mrk_1")], raw[c("id_2", "id_1"), c("mrk_2", "mrk_1")])
    expect_equal(bed[c("id_2", "id_1"), , drop = FALSE], raw[c("id_2", "id_1"), , drop = FALSE])
    expect_equal(bed[, c("mrk_2", "mrk_1"), drop = FALSE], raw[, c("mrk_2", "mrk_1"), drop = FALSE])
    expect_equal(bed[c("id_2", "id_1"), c("mrk_2", "mrk_1"), drop = FALSE], raw[c("id_2", "id_1"), c("mrk_2", "mrk_1"), drop = FALSE])

    expect_equal(bed[c("id_3", "id_1"), ], raw[c("id_3", "id_1"), ])
    expect_equal(bed[, c("mrk_3", "mrk_1")], raw[, c("mrk_3", "mrk_1")])
    expect_equal(bed[c("id_3", "id_1"), c("mrk_3", "mrk_1")], raw[c("id_3", "id_1"), c("mrk_3", "mrk_1")])
    expect_equal(bed[c("id_3", "id_1"), , drop = FALSE], raw[c("id_3", "id_1"), , drop = FALSE])
    expect_equal(bed[, c("mrk_3", "mrk_1"), drop = FALSE], raw[, c("mrk_3", "mrk_1"), drop = FALSE])
    expect_equal(bed[c("id_3", "id_1"), c("mrk_3", "mrk_1"), drop = FALSE], raw[c("id_3", "id_1"), c("mrk_3", "mrk_1"), drop = FALSE])

    # Do not modify indexes.
    i <- seq_len(6 * 3)
    bed[i]
    expect_equal(i, seq_len(6 * 3))

    i <- seq_len(6)
    j <- seq_len(3)
    bed[i, j]
    expect_equal(i, seq_len(6))
    expect_equal(j, seq_len(3))

})

test_that("length", {
    expect_equal(length(bed), length(raw))
})
