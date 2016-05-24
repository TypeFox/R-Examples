library(switchr)
removeLib("IRangesDev")
man = BiocSVNManifest()

##change this to where you want it go go
dir1 = "~/localcheckout1"
if(!file.exists(dir1))
    dir.create(dir1)
rep = lazyRepo("IRanges", man, dir = dir1)

switchTo("IRangesDev")
res = install_packages("IRanges", rep)
library("IRanges")
sessionInfo()
switchTo("IRangesDev2")
instdir = file.path(dir1, "BiocGenerics", "inst")
if(!file.exists(instdir))
    dir.create(instdir)
cat("super important dev work", file = file.path(dir1, "BiocGenerics", "inst", "myfile"))


bcdesc = file.path(dir1, "BiocGenerics", "DESCRIPTION")
desc = readLines(bcdesc)
desc = gsub("Version: .*", "Version: 1.9-9", desc)
cat(desc, sep="\n",file = bcdesc)

irdescf = file.path(dir1, "IRanges", "DESCRIPTION")
irdesc = readLines(irdescf)
irdesc = gsub("Version: .*", "Version: 4.9-9", irdesc)
irdesc = gsub("BiocGenerics \\(>= 0.13.6\\)", "BiocGenerics \\(>= 1.9.9\\)", irdesc)
cat(irdesc, sep="\n",file = irdescf)

switchTo("IRangesDev2")
rep = lazyRepo("IRanges", man, dir = dir1, force_refresh=TRUE)

install_packages("IRanges", rep)
.Internal(lazyLoadDBflush('/home/beckerg4/.switchr/IRangesDev/IRanges/R/IRanges.rdb'))
.Internal(lazyLoadDBflush('/home/beckerg4/.switchr/IRangesDev/IRanges/help/IRanges.rdb'))
