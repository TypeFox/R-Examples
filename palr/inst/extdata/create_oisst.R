# #https://github.com/AustralianAntarcticDivision/raadtools
#   library(raadtools)
#   readsst(latest = TRUE)
#   f <- tail(sstfiles(), 1)
#   oisst <- readsst(f$date, xylim = extent(140, 180, -65, -30))
#   metadata(oisst) <- list(ncdump = local({suppressWarnings(r <- raster(f$fullname)); capture.output(print(r))[-1]}), 
#  filename = gsub("data/", "", f$file))
#   save(oisst, file = "inst/extdata/oisst.rda", compress = "bzip2")
