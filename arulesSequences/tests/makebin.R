
### ceeboo 2008

library(arulesSequences)

data(zaki)

arulesSequences:::makebin(zaki, "zaki")
arulesSequences:::write_cspade(zaki, "zaki.asc")

exe <- "bin"
if (.Platform$r_arch != "")
    exe <- file.path(exe, .Platform$r_arch)
exe <- system.file(exe, package = "arulesSequences")
system2(file.path(exe, "makebin"), args = c("zaki.asc makebin.data"))

stopifnot(!system("cmp zaki.data makebin.data"))

system2(file.path(exe, "getconf"), args = c("-i makebin -o makebin"))

stopifnot(!system("cmp zaki.conf makebin.conf"))

###
