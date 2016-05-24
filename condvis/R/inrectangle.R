`%inrectangle%` <-
# helper function for checking what screen was clicked on
function (point, rectangle)
{
    # assuming (x, y) and (xleft, xright, ybottom, ytop)
    check1 <- point[1] >= rectangle[1]
    check2 <- point[1] < rectangle[2]
    check3 <- point[2] >= rectangle[3]
    check4 <- point[2] < rectangle[4]

    check1 && check2 && check3 && check4
}
