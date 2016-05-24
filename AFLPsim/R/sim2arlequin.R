sim2arlequin <-
function (x,filename) {
x$S <- NULL
x$SelMarkers <- NULL
 x <- x[!is.na(x)]
if (class(x[[1]]) == "list") {
  x <- unlist(x, recursive = F)
  x <- lapply(x, function(x) x[, -1])
}
N <- length(x)
N2 <- N 
raw <- do.call("rbind", x)
ind <- rownames(raw)
popnames <- names(x)
n <- numeric(N2)
for (i in 1:N2) {
  n[i] <- nrow(x[[i]])
}
pop <- rep(popnames[1:N2], n[1:N2])
col <- c("ind", "pop", colnames(raw))
mat <- as.data.frame(cbind(ind, pop, raw))
colnames(mat) <- col
manb <- dim(mat)[2] - 2
innb <- dim(mat)[1]
mat <- mat[1:(innb), ]
mat <- mat[order(mat[, 2]), ]
pops <- as.vector(mat[, 2])
inds <- as.vector(mat[, 1])
matm <- as.matrix(mat[, (3:(manb + 2))])
popsizes <- table(pops)
npop <- length(popsizes)
popnames <- vector(mode = "character", npop)
n <- 0
cat("[Profile]", "\n", "Title=", "\"", "AFLP data", "\"", 
    "\n", "NbSamples=", file = filename, sep = "")
cat(npop, "\n", "DataType=RFLP", "\n", "GenotypicData=0", 
    "\n", "LocusSeparator=NONE", "\n", "MissingData='?'", 
    "\n", "\n", file = filename, sep = "", append = TRUE)
cat("[Data]", "\n", "[[Samples]]", "\n", "\n", file = filename, 
    append = TRUE)
for (i in 1:npop) {
  cat("SampleName=", "\"", pops[n + 1], "\"", "\n", file = filename, 
      sep = "", append = TRUE)
  cat("SampleSize=", popsizes[i], "\n", file = filename, 
      append = TRUE)
  cat("SampleData= {", "\n", file = filename, append = TRUE)
  for (j in 1:popsizes[i]) {
    cat(inds[n + j], "\t", "1", "\t", matm[(n + j), ], 
        "\n", file = filename, append = TRUE)
  }
  cat("}", "\n", "\n", file = filename, append = TRUE)
  popnames[i] <- pops[n + 1]
  n <- n + popsizes[i]
}
cat("[[Structure]]", "\n", "\n", "StructureName = ", "\"", 
    "one group", "\"", "\n", "NbGroups = 1", "\n", file = filename, 
    sep = "", append = TRUE)
cat("\n", "Group ={", "\n", file = filename, sep = "", append = TRUE)
for (i in 1:npop) {
  cat("\"", popnames[i], "\"", "\n", file = filename, append = TRUE, 
      sep = "")
}
cat("}", "\n", file = filename, append = TRUE)
}
