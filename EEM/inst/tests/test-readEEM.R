context("test readEEM")

test_that("files read are labelled as EEM class", {
    folder <- "data_format"
    files <- c("hitachi_F-7000.txt",
               "hitachi_F-7000_Japanese.TXT",
               "JASCO_FP-8500.csv", 
               "shimadzu_RF-6000.txt",
               "horiba_aqualog.dat"
               )
    data <- readEEM(file.path(folder,files))
    expect_is(data, "EEM")
})
