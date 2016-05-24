simMixMerlin <- function (x, aa, afreq, options=NULL, seed=12345, generate = FALSE) 
{
  if(!is.null(options)) {
    merlinopt <- paste("--sim --save -r",seed,options)
  } else {
    merlinopt <- paste("--sim --save -r",seed)
  }
  merlin(x, model = F, cleanup = F, generate.files = generate, 
         options = merlinopt)
  sim.ped <- scan(file = "merlin-replicate.ped", what = c(1, 1, 1, 1, 1, rep(c("", 1)), x$nMark * 2), comment.char = "e", 
                  quiet = TRUE)
  sim.ped <- sub("/", " ", sim.ped)
  sim.ped <- matrix(sim.ped, byrow = TRUE, nrow = x$nInd)
  swap = c(1:5, ncol(sim.ped), c(6:(ncol(sim.ped) - 1)))
  sim.ped = sim.ped[, swap]
  sim.ped = matrix(as.numeric(sim.ped), nrow = x$nInd)
  write.table(sim.ped, file = "_sim.ped", quote = FALSE, col.names = FALSE, row.names = FALSE)
  y = linkdat("_sim.ped", freq = "_merlin.freq", map = "_merlin.map", dat = "_merlin.dat", verbose = FALSE)
  mixlist = lapply(y$markerdata, function(m) sort(unique(m[m > 0])))
  for (i in 1:length(mixlist)) mixlist[[i]] = aa[[i]][mixlist[[i]]]
  id = x$available[1]
  for (i in 1:x$nMark) {
    g = y$markerdata[[i]][id, 1:2]
    g = aa[[i]][g]
    y = modifyMarker(y, i, id, g, alleles = aa[[i]], afreq = afreq[[i]])
  }
  list(y = y, mixlist = mixlist)
}

