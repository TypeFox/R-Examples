f.coefnames <- function(coefnames){
##
## FIND NAMES/POSITIONS OF RELEVANT PARAMETERS. NOTE: \\< INSISTS ON START OF WORD, SO THAT, FOR INSTANCE, "cm1" ISN'T PICKED UP BY "m"
##
.mf <- grep("\\<mf\\d+\\>", coefnames, value = T)
.c <- grep("\\<c\\d+\\>", coefnames, value = T)
.cm <- grep("\\<cm\\d+\\>", coefnames, value = T)
.cf <- grep("\\<cf\\d+\\>", coefnames, value = T)
.cdd <- grep("\\<cdd\\d+\\>", coefnames, value = T)
.m <- grep("\\<m\\d+\\>", coefnames, value = T)
.mdd <- grep("\\<mdd\\d+\\>", coefnames, value = T)
.poo <- grep("\\<cm_cf\\d+\\>", coefnames, value = T)
#
.ut <- list(haplo.freq = .mf, child.s = .c, child.d = .cdd, child.poo.m = .cm, child.poo.f = .cf, maternal.s = .m, maternal.d = .mdd, poo = .poo)
#
return(.ut)
}