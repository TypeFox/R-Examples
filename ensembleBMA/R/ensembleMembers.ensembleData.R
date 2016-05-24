`ensembleMembers.ensembleData` <-
function (x) 
{ 
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 k <- attr(x, "ensembleSize")
 (dimnames(x)[[2]])[1:k]
}

