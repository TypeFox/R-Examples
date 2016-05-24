order.families <-
  function(famid, patid, fid, mid, sex = NA) {
    if (!all(is.na(sex))){
      if (! all(sex[!is.na(sex)] %in% c(0, 1, 2))) {
        stop(paste("Error in order.families: ", 
                   "sex should be 1(male), 2(female) or 0(unkown).",
                   sep = ""))
      }
      if (!all(sex[patid %in% unique(fid[fid != 0])] == 1)) {
        stop("Error in order.families: all fathers have to be sex=1")
      }
      if (!all(sex[patid %in% unique(mid[mid != 0])] == 2)) {
        stop("Error in order.families: all mothers have to be sex=2")
      }
    }
    if (any(duplicated(patid))) { 
      stop("Error in order.families: All patid values must be unique.")
    }
    if (any(patid == 0)) { 
      stop("Error in order.families: All patid have to be not equal to 0.") 
    }
    n <- length(patid)
    # check generation
    fidi <- match(fid, patid, nomatch = 0)
    midi <- match(mid, patid, nomatch = 0)
    parents <- ((fidi== 0) & (midi == 0))
    generation <- rep(0, n)
    for(i in 1:n) {
      child <- match(
        mid, 
        patid[parents], nomatch = 0) + 
        match(fid, patid[parents], nomatch = 0)
      if (!all(child == 0)) {
        if (i == n) {
          stop(paste("Error in order.families: Contradiction in pedigrees.", 
                     sep = ""))
        }
        parents <- (child > 0)
        generation[parents] <- i
      }
    }
    if (any(generation > 1)) {
      stop(paste("Error in order.families: Only nuclear families are allowed.", 
                 sep = ""))
    }
    # order families
    parent <- as.numeric(patid %in% unique(fid[fid != 0]))
    parent <- parent + 2 * as.numeric(patid %in% unique(mid[mid != 0]))
    if (!all(is.na(sex))) {
      sex <- ifelse(sex == 0 | is.na(sex), 3, sex)
      ord <- order(famid, generation, parent, sex, patid)
    } else {
      ord <- order(famid, generation, parent, patid)
    }
    return(as.integer(ord))
  }
