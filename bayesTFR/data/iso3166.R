iso3166 <- read.table(file='iso3166.txt', sep=';', header=TRUE, comment.char='#', na.strings="")
iso3166 <- iso3166[,c(1,3,4,5)]
colnames(iso3166) <- c('name', 'charcode', 'charcode3', 'uncode')
