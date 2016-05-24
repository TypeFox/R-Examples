seqLength <- function(
##title<< Determine lengths of sequences of identical elements 
##description<< The function returns a vector indicating the length of sequences of identical 
##              elements or elements of an identical type in a vector.
    array_in          ##<< numeric: input vector
    ,funct.seq = is.na##<< string or name of a function: name of the function to use. Has to return 
                      ##   TRUE for elements belonging to the sequence and FALSE otherwise
)
##details<< The function returns a vector of the same length as the input vector indicating
##          for each element in array_in how long the sequence of similar elements this value 
##          belongs to id. Zero indicates that identical values (or types) in the neighbourhood 
##          of the value. This process helps to indicate for
##          example long gaps in timeseries.
##seealso<<
##\code{\link{seqLongest}}
{
    x            <- rle(do.call(funct.seq,args=list(array_in)))
    array_result <- rep(x$lengths, x$lengths)
    array_result[!do.call(funct.seq,list(array_in))]<-0
    ##value<< numeric vector of the same length as the input vector indicating the gap length
    ##        for each single element
    names(array_result)<- NULL
    array_result
}
