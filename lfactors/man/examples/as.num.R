require(lfactors)
# create an example
let <- lfactor(4:12,
               levels=4:12,
               labels=letters[4:12])

as.numeric(let)
#same as as.numeric(4:12)
as.integer(let)
#same as 1:9

