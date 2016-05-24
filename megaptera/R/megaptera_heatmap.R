megapteraHeatmap <- function(conn, col = "gb", sort = "phylo"){
  
  close.after <- FALSE
  if ( class(conn) == "sproj" ){
    conn <- dbConnect(PostgreSQL(), host = x$db[["host"]], 
                      port = x$db[["port"]], dbname = x$db[["dbname"]], 
                      user = x$db[["user"]], password = x$db[["password"]])
    close.after <- TRUE
  }
  x <- "SELECT * FROM taxonomy JOIN locus USING (spec) ORDER BY spec"
  x <- dbGetQuery(conn, x)
  x <- x[, grep(paste("spec|_", col, sep = ""), colnames(x))]
  colnames(x) <- gsub(paste("^_", col, sep = "|_"), "", 
                      colnames(x))
  colnames(x) <- gsub("_", ".", colnames(x))
  rownames(x) <- x$spec
  x$spec <- NULL
  x <- as.matrix(x)
  x[is.na(x)] <- 0
  if ( col == "blocks" ){
    x[grep("selected", x)] <- 1
    x[grep("excluded", x)] <- 0
    mode(x) <- "numeric"
  }
  
  if (sort == "phylo" ){
    tax <- dbReadTaxonomy(conn)
    tax <- fixNodes(ladderize(tax2tree(tax)))
    x <- x[match(tax$tip.label, rownames(x)), ]
  }
  
  not <- rowSums(x) == 0
  cat("\n", length(which(not)), "species missing")
  write(rownames(x)[not], file = "missing.txt") 
  
  ## create pdf
  pdf("aaa.pdf", paper = "a4", height = 21, width = 4)
  par(mai = c(0.3, 1, .2, 0))
#   plot(tax, show.tip.label = FALSE, no.margin = TRUE)
  if ( col == "blocks" ){
    colors <- c("darkblue", "red")
  } else {
    colors <- c("darkblue", "yellow", rep("orange", 9), 
                rep("red", max(x) - 10))
  } 

  image(t(x), xaxt = "n", yaxt = "n", col = colors)
  mtext(colnames(x), side = 3, line = 0.1, at = seq(from = 0, to = 1, length.out = ncol(x)), cex = .6, 
        adj = .5, las = 1)
  rcol <- rep("black", nrow(x))
  rcol[not] <- "red"
  mtext(gsub("_", " ", rownames(x)), 
        side = 2, line = 0.1, at = seq(from = 0, to = 1, 
                                       length.out = nrow(x)), 
        cex = .35, 
        adj = 1, las = 1, col = rcol)
  dev.off()
  system("open aaa.pdf")
  system("open notOnGenBank.txt")
}