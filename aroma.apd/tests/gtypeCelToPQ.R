library("aroma.apd")

# Scan for CEL files
files <- list.files(pattern="[.](cel|CEL)$")

# Convert each to RAW file
for (file in files) {
  rawFile <- gsub("[.][^.]*$", ".raw", file)
  file.remove(rawFile)
  cel <- gtypeCelToPQ(file, verbose=TRUE)
  write.table(file=rawFile, cel, sep="\t", quote=FALSE)
}
