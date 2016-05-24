iso3166 <- utils::read.table(file='iso3166.txt', sep=';', header=TRUE, comment.char='#', na.strings="")
iso3166 <- iso3166[,c(1,3,4,5)]
colnames(iso3166) <- c('name', 'charcode', 'charcode3', 'uncode')
iso3166 <- cbind(iso3166, is.country=TRUE)

# testthat does not seem to source the file below
# source('iso3166ud.R')
# therefore this is a copy from iso3166ud.R
iso3166ud <- utils::read.table(file='iso3166ud.txt', sep=';', header=TRUE, comment.char='#', na.strings="")
colnames(iso3166ud) <- c('name', 'charcode', 'charcode3', 'uncode')
iso3166ud <- cbind(iso3166ud, is.country=FALSE)

iso3166 <- rbind(iso3166, iso3166ud)
