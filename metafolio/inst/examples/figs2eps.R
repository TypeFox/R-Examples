f <- list.files(pattern = "*.pdf")
# these have transparency, need to use Illustrator to flatten:
f <- f[-which(f == "spatial-mv.pdf")]
f <- f[-which(f == "cons-plans-n.pdf")]

for(i in 1:length(f)) {
  system(paste0("pdftops -eps ", f[i]))
}

