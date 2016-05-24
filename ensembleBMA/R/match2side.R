match2side <-
function (x, y) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 mxy <- match( x, y, nomatch = 0)

 myx <- match( y, x, nomatch = 0)

 list( xy = mxy, yx = myx)
}

