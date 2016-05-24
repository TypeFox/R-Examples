## useful functions for mapping array annotation and mutation databases to
## a reference sequence

sanitise.whitespace <- function(tt) return(gsub(" $", "", gsub("^ ", "", gsub("  ", " ", tt))))

parse.snps <- function(seqs, seqnames = NULL) {
  tmp <- as.data.frame(t(sapply(seqs, function (seq) {
    remat <- regexpr("\\[[ACGT]\\/[ACGT]\\]", toupper(seq)) # =RegExp MATch
    if (remat > 0 && attr(remat, "match.length") == 5) {
      return(c(seq5p = substr(seq, 1, remat - 1), # "" if remat==1
                  poly = substr(seq, remat, remat + 4),
                  alleleA = substr(seq, remat+1, remat + 1),
                  alleleB = substr(seq, remat+3, remat + 3),
                  seq3p = substr(seq, remat + attr(remat, "match.length"), nchar(seq)) # "" if remat+match.length > nchar(seq)
                  ))
    } else {
      return(c(seq5p = NA, poly = NA, alleleA = NA, alleleB = NA, seq3p = NA))
    }})), stringsAsFactors = FALSE)
  tmp$len5p = nchar(tmp$seq5p)
  tmp$len3p = nchar(tmp$seq3p)
  rownames(tmp) <- seqnames
  return(tmp)
}

## fastaname=tempfile() does not work with remote blat because trying to use rtp:blat//tmp//....
run.blat <- function(polys, fastaname,
                     system.pre = c("ssh rtp mkdir -p blat", "scp FASTANAME rtp:blat"),
                     system.blat = "ssh rtp blat -mask=lower -fastMap -noHead blat/dblist blat/FASTANAME blat/FASTANAME.psl",
                     system.post = "scp rtp:blat/FASTANAME.psl .",
                     dry.run = FALSE) {
  polys <- as.data.frame(polys)
  if (!dry.run) {
    if ("query" %in% names(polys)) {
      queries <- polys$query
    } else if (all(c("seq5p", "alleleA", "seq3p") %in% names(polys))) {
      queries <- paste(polys$seq5p, polys$alleleA, polys$seq3p, sep = "")
    } else {
      stop("polys must have either 'query' variable, or 'seq5p', 'alleleA' and 'seq3p' variables")
    }
    sink(fastaname)
    for (idx in 1:nrow(polys)) {
      cat(">", rownames(polys)[idx], "\n", sep = "")
      cat(queries[idx], "\n", sep = "")
    }
    sink()
    for (runme in gsub("FASTANAME", fastaname, c(system.pre, system.blat, system.post))) system(runme)
  }
  return(read.table(paste(fastaname, ".psl", sep = ""),
                    header = FALSE, 
                    col.names = c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert",
                      "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd",
                      "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"),
                    as.is = TRUE))
  ## BLAT psl header names from genome.ucsc.edu/FAQ/FAQformat#format2
}

## + strand matches, qStart and qStarts are both relative to 0 position of query sequence
## - strand matches, q[Start|End] are relative to 0 position of query sequence
##                   but qStarts is relative to 0 position of reverse complement of query sequence

remap.q2t <- function(qpos, blatres) {
  ## remap multiple positions in a *single* query (q) sequence to positions in the target (t) sequence
  blatres <- as.data.frame(blatres)
  stopifnot(nrow(blatres) == 1)
  if (blatres$blockCount == 1) {
    ## alignment is 1 block so use q[Start|End] and t[Start|End] coordinates
    if (blatres$strand == "+") {
      return(ifelse(qpos >= blatres$qStart & qpos < blatres$qEnd, qpos - blatres$qStart + blatres$tStart, NA))
    } else if (blatres$strand == "-") {
      return(ifelse(qpos >= blatres$qStart & qpos < blatres$qEnd, blatres$qEnd - 1 - qpos + blatres$tStart, NA))
    } else {
      return(rep(NA, length(qpos)))
    }
  } else {
    ## alignment is >1 blocks
    blockSizes <- as.integer(unlist(strsplit(blatres$blockSizes, ",")))
    tStarts <- as.integer(unlist(strsplit(blatres$tStarts, ",")))
    if (blatres$strand == "+") {
      qStarts <- as.integer(unlist(strsplit(blatres$qStarts, ","))) # wrong? to have blatres$qStart + 
      ## (x - qStarts) is 0-offset position within block
      return(sapply(qpos, function(x) {
        blockIdx <- which(x >= qStarts & x < qStarts + blockSizes)
        return(if (length(blockIdx) == 1) tStarts[blockIdx] + x - qStarts[blockIdx] else NA)}))
    } else if (blatres$strand == "-") {
      qEnds <- blatres$qSize - as.integer(unlist(strsplit(blatres$qStarts, ","))) - 1
      ## (qEnds - x) is 0-offset position within block
      return(sapply(qpos, function(x) {
        blockIdx <- which(x > qEnds - blockSizes & x <= qEnds)
        return(if (length(blockIdx) == 1) tStarts[blockIdx] + qEnds[blockIdx] - x else NA)}))
    } else {
      return(rep(NA, length(qpos)))
    }
  }
}





## BAD does not deal with blockcount>1
fooo <- function(polys, blatout) {
  hitcount <- table(blatout$Q.name)
  hits <- as.vector(hitcount[match(rownames(polys), names(hitcount))])
  hits[is.na(hits)] <- 0

  stopifnot(all(blatout$strand %in% c("+", "-")))
  pos1 <- ifelse(blatout$strand == "+", blatout$tStart, blatout$T.end - 1) +
    ifelse(blatout$strand == "+", 1, -1) *
      (polys$len5p[match(blatout$Q.name, rownames(polys))] - blatout$qStart)
  tmp <- data.frame(hits = hits, pos = as.vector(ifelse(hits == 1, pos1[match(rownames(polys), blatout$Q.name)], NA)))
  rownames(tmp) <- rownames(polys)
  return(tmp)
}
  


subsequences <- function(pos, seq, left = 30, right = 30) { # extract base at pos and flanking sequence
  tmp <- as.data.frame(t(sapply(unique(sort(pos)), function(thispos) {
    if (is.na(thispos) || thispos < 1 || thispos > nchar(seq)) return(c(pos = thispos, base = NA, leftLen = NA, rightLen = NA, subsequence = NA))
    thisleft <- if (thispos - left < 1) thispos - 1 else left 
    thisright <- if (thispos + right > nchar(seq)) nchar(seq) - thispos else right
    return(c(pos = thispos, base = substr(seq, thispos, thispos), leftLen = thisleft, rightLen = thisright, subsequence = substr(seq, thispos - thisleft, thispos + thisright))) # R's substr uses 1-offset, +0 end
  })), stringsAsFactors = FALSE)
  for (colname in c("pos", "leftLen", "rightLen")) tmp[[colname]] <- as.integer(tmp[[colname]])
  return(tmp)
}
