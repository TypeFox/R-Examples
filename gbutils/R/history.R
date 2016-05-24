                                                         # moved from psistudent on 2014-02-20
raw_history <- function()
{
    file1 <- tempfile("Rrawhist")
    savehistory(file1)
    rawhist <- readLines(file1)
    unlink(file1)

    rawhist
}
