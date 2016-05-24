iso3166ud <- utils::read.table(file='iso3166ud.txt', sep=';', header=TRUE, comment.char='#', na.strings="")
colnames(iso3166ud) <- c('name', 'charcode', 'charcode3', 'uncode')
iso3166ud <- cbind(iso3166ud, is.country=FALSE)