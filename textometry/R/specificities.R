#* Copyright © - 2008-2011 ANR Textométrie - http://textometrie.ens-lyon.fr
#* Copyright © - 2012-2014 ENS de Lyon - http://textometrie.ens-lyon.fr
#*
#* This file is part of the TXM platform - http://sourceforge.net/projects/txm.
#*
#* The TXM platform is free software: you can redistribute it and/or modify
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* The TXM platform is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#* General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with the TXM platform.  If not, see <http://www.gnu.org/licenses/>.

## Sylvain Loiseau <sloiseau@u-paris10.fr>
## Lise Vaudor <lise.vaudor@ens-lyon.fr>
## Matthieu Decorde <matthieu.decorde@ens-lyon.fr>
## Serge Heiden <slh@ens-lyon.fr>

# export TOC
# phyper_right							Calculate specificity probabilities (call dhyper)
# specificities                         Calculate Lexical Specificity Score (call specificities.probabilities)
# specificities.banals					*PROTOTYPE* specificity lines between max and min banality score thresholds
# specificities.distribution.plot       Display specificities probability distribution (call dhyper and specificities.probabilities.vector)
# specificities.lexicon                 Specificities association score with two frequency list (call specificities)
# specificities.lexicon.new             *OBSOLETE* (compatibility holder, see 'specificities.lexicon')
# specificities.probabilities           Calculate specificity probabilities on lexical table (call specificities.probabilities.vector)
# specificities.probabilities.vector	Calculate specificity probabilities on vector (call phyper and phyper_right)
# specificities.top						*PROTOTYPE* top n max or min specificity lines

`specificities.distribution.plot` <- function(x, F, t, T) {
  
  # x: observed number of A words
  # F: total number of A
  # t: size of part
  # T: size of corpus
  
  # return : (px, mode, pfsum)
  
  # example: specificities.distribution.plot(11, 296, 1084, 61449)
  
  f=0:F
  
  mode=((F+1)*(t+1))/(T+2)

  xRightLimit=(F/4)
  if (xRightLimit < mode) xRightLimit = (F/4) + mode
  rightDelta=xRightLimit/10
  posType=4
  
  pf=dhyper(f, F, T-F, t)
  pfsum=specificities.probabilities.vector(f, rep(F, F+1), rep(T, F+1), t)
  
  px <- pf[x+1]
  pfsumx <- pfsum[x+1]

  aY=max(pf)/3
  
  # window width around mode
  thr <- 1.0e-20
  i_min=0
  i_max=0
  for (i in mode:0) {
    if (pf[i+1] < thr) {
      i_min = i
      break
    }
  }
  for (i in mode:F) {
    if (pf[i+1] < thr) {
      i_max = i
      break
    }
  }
  w <- abs(i_max - i_min + 1)
  
  # x coordinates to display
  x1 <- max(min(mode, x) - w, 0)
  x2 <- min(max(mode, x) + w, F)
  
  # arrow and text coordinates
  if (aY < px) {
    obsY = aY*1.3
  } else if (x < mode) {
    obsY = aY*1.3
  } else {
    obsY = aY*0.7
  }
  
  # plot(f, pfsum, xlab="f'", ylab="P", type="l", col="blue") c(x1, x2, 0, max(px))
  plot(f, pf, ylab="Pspecif", type="l", col="green", xlim=c(x1, x2))
  Axis(at=c(mode, x), side=1, labels=c(sprintf("%d", round(mode)), sprintf("%d", x)))

  arrows(x+rightDelta, aY, mode, 0, length=0.1)   # mode arrow
  arrows(x+rightDelta, obsY, x, px, length=0.1) # obs arrow
  
  smode <- sprintf("mode = %d", round(mode))
  spf <- sprintf("p(f=%d) = %s\np(f >= %d) = %s", x, format.default(px, scientific = TRUE, digits = 4), x, format.default(pfsumx, scientific = TRUE, digits = 4))
  
  text(x+rightDelta, aY, smode, pos=posType, xpd=TRUE)
  text(x+rightDelta, obsY, spf, pos=posType, xpd=TRUE)
 
  title(sprintf("Distribution de probabilit\u00E9 de la sp\u00E9cificit\u00E9\nde param\u00E8tres F=%d, t=%d, T=%d,\nfr\u00E9quence observ\u00E9e = %d", F, t, T, x))
  
  # for cumulative distribution plot
  # s2 <- sprintf("integral P(f' >= f)")
  # legend("topright", c(str1, s2),lty = 1, col=c("green", "blue"), inset = .02, cex=0.9)
  # legend("topright", c(str1), lty = 1, col=c("green"), inset = .02, cex=0.9)
  
  return(list(px=px, mode=mode, pfsum=pfsum))
}

println <- function(prefix, obj) {
  cat(prefix)
  print(obj)
}

phyper_right=function(v_white_drawn, v_num_white, v_num_black, num_drawn){

# println("v_white_drawn: ", length(v_white_drawn))
# println("v_num_white: ", length(v_num_white))
# println("v_num_black: ", length(v_num_black))
# println("num_drawn: ", num_drawn)
  
  v_s2=rep(NA,length(v_white_drawn))

  for (j in 1:length(v_white_drawn)){
    a=v_white_drawn[j]
    b=v_num_white[j]
    c=v_num_black[j]
    d=num_drawn
    a_tmp=a#+1
    s1=dhyper(a_tmp, b, c, d)
    s_tmp=dhyper(a_tmp+1, b, c, d)
    s2=s1+s_tmp
    a_tmp=a_tmp+1
    while(log(s2)!=log(s1)){
      s1=s2
      a_tmp=a_tmp+1
      s_tmp=dhyper(a_tmp, b, c, d)
      s2=s1+s_tmp
    }
    v_s2[j]=s2
  }
  return(v_s2)
}

# old function format for retro-compatibility
`specificites` <-
  function(lexicaltable, types=NULL, parts=NULL) {
    return(specificities(lexicaltable, types, parts));
  }

`specificities` <-
  function(lexicaltable, types=NULL, parts=NULL) {
    spe <- specificities.probabilities(lexicaltable, types, parts);
    spelog <- matrix(0, nrow=nrow(spe), ncol=ncol(spe));
    spelog[spe < 0] <- log10(-spe[spe < 0]);
    spelog[spe > 0] <- abs(log10(spe[spe >0]));
    spelog[spe == 10] <- +Inf;
    spelog[spe == -10] <- -Inf;
    #spelog[is.infinite(spe)] <- 0; # cannot be infinite
    spelog <- round(spelog, digits=4);
    rownames(spelog) <- rownames(spe);
    colnames(spelog) <- colnames(spe);
    #class(spelog) <- "specificities";
    #attr(spelog, "l.t") <- spe;
    return(spelog);
  }

`specificities.probabilities` <-
  function(lexicaltable, types=NULL, parts=NULL) {
    
    # if (!is.numeric(lexicaltable)) stop("The lexical table must contain numeric values.");
    
    rowMargin <- rowSums(lexicaltable); # or "F" (the total frequency of all the types).
    colMargin <- colSums(lexicaltable); # or "t" (the size of the parts).
    F <- sum(lexicaltable);             # or "T" (the grand total, number of tokens in the corpus).
    
    if (! is.null(types)) {      # Filter on tokens to be considered.
      if(is.character(types)) {  # convert the name of types given with "types" into row index numbers.
        if (is.null(rownames(lexicaltable))) stop("The lexical table has no row names and the \"types\" argument is a character vector.");
        if (! all(types %in% rownames(lexicaltable))) stop(paste(
          "Some requested types are not known in the lexical table: ",
          paste(types[! (types %in% rownames(lexicaltable))], collapse=" "))
        ); 
      } else {
        if (any(types < 1)) stop("The row index must be greater than 0.");
        if (max(types) > nrow(lexicaltable)) stop("Row index must be smaller than the number of rows.");
      }
      lexicaltable <- lexicaltable[types, , drop = FALSE];
      rowMargin <- rowMargin[types];
    }
    
    if (! is.null(parts)) {      # Filter on parts to be considered.
      if(is.character(parts)) {  # convert the name of parts given with "parts" into col index numbers.
        if (is.null(colnames(lexicaltable))) stop("The lexical table has no col names and the \"parts\" argument is a character vector.");
        if (! all(parts %in% colnames(lexicaltable))) stop(paste(
          "Some requested parts are not known in the lexical table: ",
          paste(parts[! (parts %in% colnames(lexicaltable))], collapse=" "))
        ); 
      } else {
        if (max(parts) > ncol(lexicaltable)) stop("Column index must be smaller than the number of cols.");
        if (any(parts < 1)) stop("The col index must be greater than 0.");
      }
      lexicaltable <- lexicaltable[,parts, drop=FALSE];
      colMargin <- colMargin[parts];
    }
    
    if (nrow(lexicaltable) == 0 | ncol(lexicaltable) == 0) {
      stop("The lexical table must contains at least one row and one column.");
    }
    
    specif <- matrix(0.0, nrow=nrow(lexicaltable), ncol=ncol(lexicaltable));
    
    for(i in 1:ncol(lexicaltable)) {    # We proceed the whole lexical table by column (i.e. by part).
      
      specif[,i] <- specificities.probabilities.vector(lexicaltable[,i], rowMargin, F, colMargin[i])
    }
    
    colnames(specif) <- colnames(lexicaltable);
    rownames(specif) <- rownames(lexicaltable);
    
    return(specif);
  }

`specificities.probabilities.vector` <- function(v_f, v_F, T, t) {
  
  v_whiteDrawn <- v_f;  # The frequencies observed in this part for each type.
  v_white <- v_F;       # The total frequencies in the corpus for each type.
  v_black <- T-v_white; # The total complement frequency in the corpus for each type.
  drawn <- t;           # The number of tokens in the part.
  
  v_specif <- c()
  
  v_independance    <- (v_white * drawn) / T;         # The theoretic frequency of each type.
  v_specif_negative <- v_whiteDrawn <  v_independance;  # index of observed frequencies below the theoretic frequencies.
  v_specif_positive <- v_whiteDrawn >= v_independance;  # index of observed frequencies above the theoretic frequencies.
  
  v_negativeScores <- -phyper (
    v_whiteDrawn[v_specif_negative], v_white[v_specif_negative], v_black[v_specif_negative], drawn
  );
  v_negativeScores[v_negativeScores==0] = -10
  v_specif[v_specif_negative] <- v_negativeScores
  
  v_positiveScores <- phyper_right (
    v_whiteDrawn[v_specif_positive], v_white[v_specif_positive], v_black[v_specif_positive], drawn
  );
  
  v_positiveScores[v_positiveScores==0] = 10
  v_specif[v_specif_positive] <- v_positiveScores
  
  return(v_specif)
}

`specificities.lexicon` <- function(lexicon, sublexicon) {
  lexicaltable <- lexiconsToLexicalTable(lexicon, sublexicon);
  return(specificities(lexicaltable, NULL, NULL));
}

`lexiconsToLexicalTable` <- function(lexicon, sublexicon) {
  if (! all(names(sublexicon) %in% names(lexicon))) 
    stop(paste(
      sum(! (names(sublexicon) %in% names(lexicon))),
      "types of the sublexicon not found in the lexicon: ",
    )
    ); 
  
  sub <- numeric(length(lexicon));
  names(sub) <- names(lexicon);
  sub[names(sublexicon)] <- sublexicon;
  
  complementary.lexicon <- c(lexicon- sub);
  if (any(complementary.lexicon < 0)) 
    stop("type cannot be more frequent in the sublexicon than in the lexicon");
  
  lexicaltable <- matrix(c(sub, complementary.lexicon), ncol=2);
  rownames(lexicaltable) <- names(lexicon);
  colnames(lexicaltable) <- c("sublexicon", "complementary");
  return(lexicaltable)
}

`specificities.lexicon.new` <- function(lexicon, sublexicon) {
  return(specificities.lexicon(lexicon, sublexicon))
}

`specificities.top` <- function(specif, N=10) {
  top <- c()
  cols <- colnames(specif)
  
  for (i in 1:length(cols)) {
    sorted <- sort(specif[, cols[i]], decreasing=TRUE)
    top <- union(top, names(sorted[1:N]))
    len <- length(sorted)
    top <- union(top, names(sorted[len-N: len]))
  }
  
  return(top)
}

`specificities.banals` <- function(specif, score=1.0) {
  banals <- rownames(specif)
  cols <- colnames(specif)
  for (i in 1:length(cols)) {
    col <- specif[, i]
    col <- col[col > score | col < -score]
    banals <- setdiff(banals, names(col))
  }
  
  return(banals)
}

#print.specificities(x, line=20, part=1, form=NULL, ...) {
#  if (all(is.null(line, part))) {
#    stop("either a line or a part must be specified");
#  }
#  if (all(!is.null(line, part))) {
#    stop("only a line or a part must be specified");
#  }
#}
