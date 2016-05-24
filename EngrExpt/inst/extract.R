dsn <- c("absorb", "adhesion2", "adhesion", "alum", "applicat",
         "assay", "bacteria", "bath", "battery", "break", "bright",
         "calcium", "caliper", "ccthickn", "cement", "cheese",
         "chemreac", "computer", "cure", "curl", "defoam", "deink2",
         "deink", "dhaze", "diagnostic", "diameter", "drought",
         "drums", "dry", "epoxy", "exposure", "fbuild", "fill",
         "fillweight", "fish2", "fish", "fluoride", "gloss", "labcomp",
         "lwsw", "lw", "moisture", "mw", "odor", "oven", "phmeas",
         "ph", "pigment", "protein", "purity", "railcar2", "railcar3",
         "railcar", "ratings2", "ratings", "reflect", "safety", "sales",
         "sarea", "separate", "soap", "stab", "stretch", "surfarea",
         "tablets", "temprate", "tennis", "thinfilm", "timetemp",
         "tpaste", "urine", "uvcoatin", "uvoven", "viscosity",
         "vitamin", "wash", "water", "webtraff", "webvisit",
         "weight", "whitearea", "yellow", "yield")
zfile <- system.file("TXT.zip", package = "EngrExpt")
pkgbase <- "~/src/R-forge/EngrExpt/pkg"
factorize <- function(x) {
    ndistinct <- length(unique(x))
    if (is.integer(x) && ndistinct < 6 && all(x %in% 0:ndistinct))
        return(factor(x, labels = LETTERS[seq_len(ndistinct)]))
    return(x)
}
for (nm in dsn) {
    if (nm %in% c("bath", "exposure", "fill", "protein", "stretch")) {
        tmp <- read.table(unz(zfile,
                              paste("TXT/", nm, ".txt", sep = '')),
                          header = TRUE)
    } else {
        tmp <- read.csv(unz(zfile,
                              paste("TXT/", nm, ".txt", sep = '')))
    }
    tmp <- do.call(data.frame, lapply(tmp, factorize))
    names(tmp) <- tolower(names(tmp))
    assign(nm, tmp)
    save(list = nm, file = file.path(pkgbase, "data", paste(nm, ".rda", sep = '')))
#    prompt(name = nm, filename = file.path(pkgbase, "man", paste(nm, ".Rd", sep = '')))
}
rm(tmp)
ls.str()
