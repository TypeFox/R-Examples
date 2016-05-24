## For a dna_seg, auto-annotates:
## Takes gene name if possible, nothing or locus_tag for the rest.
## For operons or sequences of genes, creates spanning annotations
## Gene names have to be consecutive and end with a number or a capital letter
auto_annotate <- function(dna_seg, locus_tag_pattern=NULL, names=dna_seg$gene,
                          keep_genes_only=TRUE, ...){
  nr <- nrow(dna_seg)
  if (is.null(names) || length(names) != nr)
    stop("names should be of the same length as dna_seg")
  idx <- rep(TRUE, nr)
  x1 <- middle(dna_seg)
  x2 <- rep(NA, nr)
  last_num <- numeric(nr)
  last_char_idx <- numeric(nr)
  prefix <- character(nr)
  for (i in 1:nr){
    ## If gene is not named
    if (names[i] == "-" || names[i] == ""){
      if (keep_genes_only){
        idx[i] <- FALSE
      } else if (!is.null(locus_tag_pattern)
                 && is.character(locus_tag_pattern)) {
        names[i] <- gsub(paste("^", locus_tag_pattern, "0+(.*)", sep=""),
                         "\\1", dna_seg$name[i], ignore.case=TRUE)
      } else {
        names[i] <- dna_seg$name[i]
      }
    } else {
      ## Else, which "type" is the name
      type <- NA
      if (length(grep("[A-Z]$", names[i])) > 0){
        last_char_idx[i] <- which(LETTERS == gsub(".*(.$)", "\\1", names[i]))
        prefix[i] <- gsub("^(.*).$", "\\1", names[i])
      } else if (length(grep("\\d+$", names[i])) > 0){
        ln <- gsub("^\\D*(\\d+$)", "\\1", names[i])
        last_num[i] <- if (ln != names[i]) ln else 0
        prefix[i] <- gsub("^(\\D+)\\d+$", "\\1", names[i])
      }
    }
  }
  #browser()
  ## Make sure last_* are numeric
  last_num <- as.numeric(last_num)
  ## Letter series
  series_idx <- numeric(0)
  series_dir <- 0
  series_type <- "" # num or char
  for (i in 1:(nr+1)){
    break_me <- FALSE
    ## Has pref + last_num/last_char?
    if (i <= nr && nchar(prefix[i]) > 0 &&
          (last_char_idx[i] > 0 || last_num[i] > 0)){
      valid <- TRUE
      curr_type <- if (last_num[i] > 0) "num" else "char"
      ## Is series started
      if (length(series_idx) > 0){
        ## Is type + prefix compatible?
        if (prefix[i] == prefix[series_idx[1]] && curr_type == series_type){
          ## Has series direction?
          if (series_dir != 0){
            ## Can I increment?
            if ((series_type == "num" &&
                 last_num[i-1] + series_dir == last_num[i]) ||
                (series_type == "char" &&
                 last_char_idx[i-1] + series_dir == last_char_idx[i])){
              ## Increment (do nothing)
            } else {
              ## Else break
              break_me <- TRUE
            }
          } else {
            ## No direction. Can I give one?
            dif <- ifelse(series_type == "num", last_num[i] - last_num[i-1],
                          last_char_idx[i] - last_char_idx[i-1])
            if (abs(dif) == 1){
              series_dir <- dif
            } else {
              ## Else break
              break_me <- TRUE
            }
          }
        } else {
          ## Type + prefix not compatible
          break_me <- TRUE
        }
      } else {
        ## Series not started, start it
        series_type <- if (last_num[i] > 0) "num" else "char"
      }
    } else {
      ## No pref + last_num/last_char
      break_me <- TRUE
      valid <- FALSE
    }
    ## Break or increment series
    if (break_me || i > nr){
      if (length(series_idx) > 1){
        ## Record series
        idx[series_idx[2:length(series_idx)]] <- FALSE
        x1[series_idx[1]] <- dna_seg$start[series_idx[1]]
        x2[series_idx[1]] <- dna_seg$end[series_idx[length(series_idx)]]
        suffices <- if(series_type == "num") {
          last_num[series_idx]
        } else {
          LETTERS[last_char_idx[series_idx]]
        }
        names[series_idx[1]] <- paste(prefix[series_idx[1]], suffices[1], "-",
                                      suffices[length(suffices)], sep="")
      } else {
        ## No series to record, thus do nothing
      } 
      ## Reset series
      series_dir <- 0
      ## Start a new one at once if valid
      if (valid){
        series_idx <- i 
        series_type <- curr_type
      } else {
        series_idx <- numeric(0)
        series_type <- ""
      }
    } else {
      ## Else increment
      series_idx <- c(series_idx, i)
    }
  }
  annots <- annotation(x1=x1, x2=x2, text=names, ...)
  annots[idx,]
}
