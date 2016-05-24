`contiguous.zoo` <- function(x)
{
z.rle <- rle(!is.na(rowSums(coredata(x))))
# row indexes
ends <- cumsum(z.rle$lengths)
starts <- ends - z.rle$lengths + 1
indexes <- with(z.rle, data.frame(starts, ends, lengths, values))
indexes.sort <- indexes[order(-indexes$lengths), ]
indexes.sort[indexes.sort$values, ]
}