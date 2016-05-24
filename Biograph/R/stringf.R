stringf <-
function (string)   # D:\Data\DEA HRS\HRS data and analysis\Analysis\R MR\
#  Converts a character string into a vector of characters (str_char)
# Example stringf("test") => "t" "e" "s" "t"
# SEE ALSO:  substring("statistics",1:nchar("statistics"),1:nchar("statistics"))

{   numchar <- nchar(string)
    str_char <- vector(mode="character",length=numchar)
    for (k in 1:numchar) str_char[k] <- substr(string,k,k)
    return  (str_char)
}
