read.documents <-
function (filename = "mult.dat") 
{
    one <- scan(filename, what = "", sep = "\n")
    two <- chartr(":", " ", one)
    three <- strsplit(two, " ", fixed = TRUE)
    lapply(three, function(x) matrix(as.integer(x[-1]), nrow = 2))
}
