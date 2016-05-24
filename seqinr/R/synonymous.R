syncodons <- function(codons, numcode = 1) 
{
  #
  # Check argument:
  #
  if( any(nchar(codons) != 3)) stop("vector of three character string elements expected")

  #
  # Force to lower case:
  #
  alph <- unlist(lapply(codons, s2c))
  if(any(alph %in% LETTERS)) codons <- tolower(codons)

  #
  # Get the genetic code map:
  #
  allaminos <- sapply(words(), function(x) translate(s2c(x), numcode = numcode))

  getsyn <- function(c) names(allaminos[allaminos == allaminos[[c]]])
  synonymous <- lapply(codons, getsyn)
  names(synonymous) <- codons
  return(synonymous)
}

synsequence <- function (sequence, numcode = 1, ucoweight = NULL) 
{
    tra = translate(sequence, numcode = numcode)
    cod = splitseq(sequence)
    if (is.null(ucoweight)) {
        for (a in unique(tra)) {
            pos = which(tra == a)
            if (length(pos) > 1) {
                newcod = sample(cod[pos])
                cod[pos] = newcod
            }
        }
    }
    else {
        for (a in unique(tra)) {
            pos = which(tra == a)
            urne = rep(names(ucoweight[[a]]), ucoweight[[a]] * 
                length(pos))
            if (length(urne) > 1) {
                newcod = sample(urne, length(pos))
            }
            else if (length(urne) == 1) {
                newcod = urne
            }
            else {
                print(a)
                stop("bad codon usage content")
            }
            cod[pos] = newcod
        }
    }
    newseq = s2c(c2s(cod))
    return(newseq)
}

ucoweight <- function (sequence, numcode = 1) 
{
    allaminos = s2c(c2s(SEQINR.UTIL$CODES.NCBI$CODES[numcode]))
    allcodons = splitseq(as.vector(t(cbind(rep(s2c("tcag"), each = 16), 
        rep(s2c("tcag"), each = 4), rep(s2c("tcag"), 4)))))
    syncodons = lapply(seq(21), function(a) {
        which(allaminos == unique(allaminos)[a])
    })
    usage = uco(sequence)[allcodons] #ré-ordonner selon NCBI
    weight = lapply(seq(21), function(b) {
        usage[syncodons[b][[1]]]
    })
    names(weight) = unique(allaminos)
    return(weight)
}
