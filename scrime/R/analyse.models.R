`analyse.models` <-
function (file, size.freq = TRUE, moco = c(20, 10), int.freq = TRUE,
    kmax = 10, int.level = 2, bin.names = NULL)
{
    ones <- function(model) {
        (model[begin] == 1) * abs(model[begin + 1])
    }
    twos <- function(model) {
        ifelse(model[begin] == 2, paste(model[begin + 1], model[begin +
            2]), "0 0")
    }
    threes <- function(model) {
        ifelse(model[begin] == 3, paste(model[begin + 1], model[begin +
            2], model[begin + 3]), "0 0 0")
    }
    sub.names <- function(obj, new.names) {
        nbin <- length(new.names)
        temp <- paste("", names(obj), "")
        for (i in -nbin:nbin){
            temp <- gsub(pattern = paste("", i, ""), replacement = paste("",
                paste(ifelse(i < 0, "-", ""), new.names[abs(i)], sep = ""), ""), x = temp)
         }
        end <- unlist(lapply(strsplit(temp, split = ""), length)) - 1
        temp <- substring(temp, first = 2, last = end)
        temp
    }
    erg <- read.table(file)
    n <- dim(erg)[1]
    fac <- ifelse(int.freq, 1, n)
    sizes <- ifelse(size.freq, table(erg[, 1]), table(erg[,1])/n)
    models <- erg[, 2:((int.level + 1) * kmax + 1)]
    begin <- seq(1, by = int.level + 1, length.out = kmax)
    only.ones <- apply(models, 1, ones)
    only.ones <- sort(table(as.vector(only.ones[only.ones > 0])), decreasing = TRUE)/fac
    if (!is.null(bin.names))
        names(only.ones) <- sub.names(only.ones, new.names = bin.names)
    only.twos <- apply(models, 1, twos)
    only.twos <- sort(table(as.vector(only.twos[only.twos !=
        "0 0"])), decreasing = TRUE)/fac
    if (!is.null(bin.names)) {
        names(only.twos) <- sub.names(only.twos, new.names = bin.names)
    }
    if (length(moco) > 2) {
        only.threes <- apply(models, 1, threes)
        only.threes <- sort(table(as.vector(only.threes[only.threes !=
            "0 0 0"])), decreasing = TRUE)[1:moco[3]]/fac
        if (!is.null(bin.names)) {
            names(only.threes) <- sub.names(only.threes, new.names = bin.names)
        }
    }
    else {
        only.threes <- NULL
    }
    list(size = sizes, ones = only.ones[1:moco[1]], twos = only.twos[1:moco[2]],
        threes = only.threes)
}

