PchList <- function(n)
{
    # Designed to be used in conjunction with colours to provide unambiguous
    # and rapid identification of points
    long_list <- c(                        # square, circle, triangle, diamond
        15, 16, 17, 18,                    # as solid
        0, 1, 2, 5,                        # as open
        7, 13, 9                           # as open with cross (no triangle)
    )
    if (n <= 11) {
        pch_list <- long_list[1:n]
    } else {
        # recycle the list (the colour list is of a different length)
        # long_list has a length of 11, col_list is 15
        temp_pch <- rep(long_list, ceiling(n / 11))
        pch_list <- temp_pch[1:n]
    }
    return(pch_list)
}
