# load(file.choose())
# write.table(cbind(1:281, nydf$latitude, nydf$longitude), 
#             file = "~/Desktop/coords.txt", sep = " ", row.names = FALSE, col.names = FALSE)
# write.table(cbind(1:281, 2000, nydf$population), file = "~/Desktop/pop.txt", sep = " ", col.names = FALSE, row.names = FALSE)
# write.table(cbind(1:281, floor(nydf$cases), 2000), file = "~/Desktop/cases.txt", sep = " ", col.names = FALSE, row.names = FALSE)
# 
# library(smerc)
# data(nydf)
# data(nyw)
# 
# setwd("~/Desktop")
# write.table(cbind(1:281, nydf$latitude, nydf$longitude), file = "flex_nycoords.txt", 
#             sep = " ", col.names = FALSE, row.names = FALSE)
# 
# for(i in 1:281)
# {
#   nb_i = which(nyw[i, ] == 1)
#   write(t(c(i, nb_i)), file = "flex_nyw.txt", 
#               sep = " ", ncolumns = length(nb_i) + 1, append = TRUE)
# }
# 
# write.table(cbind(1:281, floor(nydf$cases), sum(floor(nydf$cases))/sum(nydf$population)*nydf$population),             file = "flex_nycases.txt", 
#             sep = " ", col.names = FALSE, row.names = FALSE)
