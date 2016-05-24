"print.knncat" <- 
function (x, ...) 
{
#
# Print.knncat: print method for knncat objects.
#
#  Arguments: x: a knncat classifier, from knncat()
#
# Right now all this does is print the error rate.
#
yes <- sum (diag(x$misclass.mat))
tot <- sum (x$misclass.mat)
#
# Get an upper-case letter to start. Neatness counts.
#
cute <- ifelse (x$misclass.type == "train", "Training", "Test")
cat (paste (cute, " set misclass rate: ", 
    round (100 * (1 - yes/tot), 2), "%\n", sep=""))
}
