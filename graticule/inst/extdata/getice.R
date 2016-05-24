icefile <- "ftp://sidads.colorado.edu/pub/DATASETS/nsidc0051_gsfc_nasateam_seaice/final-gsfc/south/daily/2014/nt_20140320_f17_v01_s.bin"
tfile <- file.path("inst", "extdata", basename(icefile))
if (!file.exists(tfile)) download.file(icefile, tfile, mode = "wb")
