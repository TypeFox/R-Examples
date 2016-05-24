uco <- function (seq, frame = 0, index = c("eff", "freq", "rscu"), 
as.data.frame = FALSE, NA.rscu = NA) 
{
    choice <- match.arg(index)
    
    if(any(seq%in%LETTERS)){
      seq <- tolower(seq)
    }
    sequence <- splitseq(seq = seq, frame = frame, word = 3)

    if( as.data.frame == FALSE ) {
      eff <- table(factor(sequence, levels = SEQINR.UTIL$CODON.AA$CODON))
      if(choice == "eff") return(eff)
      
      freq <- eff/(floor(length(seq)/3))
      if(choice == "freq") return(freq)
      
      T <- split(freq, SEQINR.UTIL$CODON.AA$AA)
      rscu <- lapply(T, function(x) {
        return(x/((1/length(x)) * sum(x)))
      })
      names(rscu) <- NULL
      rscu <- unlist(rscu)[as.character(SEQINR.UTIL$CODON.AA$CODON)]
      is.na(rscu[!is.finite(rscu)]) <- TRUE
      rscu[is.na(rscu)] <- NA.rscu
      return(rscu)
    } else { # return all indices in a data.frame
      eff <- table(factor(sequence, levels = SEQINR.UTIL$CODON.AA$CODON))
      freq <- eff/(floor(length(seq)/3))
      T <- split(freq, SEQINR.UTIL$CODON.AA$AA)
      rscu <- lapply(T, function(x) {
        return(x/((1/length(x)) * sum(x)))
      })
      names(rscu) <- NULL
      rscu <- unlist(rscu)[as.character(SEQINR.UTIL$CODON.AA$CODON)]
      is.na(rscu[!is.finite(rscu)]) <- TRUE
      rscu[is.na(rscu)] <- NA.rscu
      df <- data.frame(SEQINR.UTIL$CODON.AA$AA, eff = eff, freq = as.vector(freq), RSCU = rscu)
      names(df) = c("AA", "codon", "eff", "freq", "RSCU")
      return(df)
    }
}


dotchart.uco <- function(x, numcode = 1, aa3 = TRUE, cex = 0.7, 
  alphabet = s2c("tcag"), pch = 21, gpch = 20, bg = par("bg"), 
  color = par("fg"), gcolor = par("fg"), lcolor = "gray", xlim, ...)
{
  if( is.null(names(x)) ) names(x) <- words( alphabet = alphabet )
#
# General sorting 
#
  x <- sort(x)
  labels <- names(x)
  stringlabel = paste(labels, sep="", collapse="")
  groups <- as.factor(translate(s2c(stringlabel), numcode =  numcode))
  gdata <- sapply(split(x, groups), sum)
#
# Now, sorting by aa order
#
  gordered <- rank(gdata)
  xidx <- numeric(64)

  for( i in seq_len(64) )
  {
    xidx[i] <- -0.01*i + gordered[groups[i]]
  }

  x <- x[order(xidx)]
  labels <- names(x)
  stringlabel = paste(labels, sep="", collapse="")
  aa <- translate(s2c(stringlabel), numcode =  numcode)
  groups <- factor(aa, levels = unique(aa))
  gdata <- sapply(split(x, groups), sum)

  if( missing(xlim) ) xlim <- c(0, max(gdata))
  if( aa3 )
  {
    levels(groups) <- aaa(levels(groups))
  }
  dotchart(x = x, labels = labels, groups = groups, gdata = gdata,
   cex = cex, pch = pch, gpch = gpch, bg = bg, color = color,
   gcolor = gcolor, lcolor = lcolor, xlim, ...)
#
# Return invisibly for further plots
#
  result <- list(0)
  result$x <- x
  result$labels <- labels
  result$groups <- groups
  result$gdata <- gdata

  ypg <- numeric( length(levels(groups)) )
  i <- 1
  for( aa in levels(groups) )
  {
    ypg[i] <- length(which(groups == aa)) + 2
    i <- i + 1
  }
  ypg <- rev(cumsum(rev(ypg))) - 1
  names(ypg) <- levels(groups)
  result$ypg <- ypg

  ypi <- numeric( length(x) )
  for( i in seq_len(length(x)) )
  {
    ypi[i] <- ypg[groups[i]]
  }
  antirank <- function(x) 
  {
    return( seq(length(x),1,by=-1 ))
  }
  ypi <- ypi - unlist(sapply(split(x, groups),antirank))
  names(ypi) <- labels
  result$ypi <- ypi

  return( invisible(result) ) 
}



