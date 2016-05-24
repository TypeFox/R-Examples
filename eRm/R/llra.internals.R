#internal functions
get_item_cats <- function(X, nitems, grp_n){
  # returns list of vectors with length max(categories) for each item;
  # 1:number categories are the first few entries and the rest is filed with zeros
  # This later corresponds to the necessary setup in LPCM where the superfluous categories must be set to 0
  its     <- rep(seq_len(nitems), each = sum(grp_n))
    ###mjm fix 2014-09-24: split() works differently with matrices (column-wise) and data.frames (row-wise), so: as.matrix() to be sure.
    ###mjm fix 2014-09-24: if NAs in raw data X, results would be NA, so: added na.rm = TRUE to the routine
  cats    <- lapply(split(as.matrix(X), its), max, na.rm = TRUE) #splits the data matrix according to items and finds the maximum category
    ###mjm fix 2014-09-24: fix for NAs, see above
  max.cat <- max(X, na.rm = TRUE) #overall maximum category
  vec.cat <- lapply(cats, function(x){ c(seq_len(x), rep(0, max.cat-x)) })
  vec.cat #the ominous list of form c(1:categories, 0, 0, 0)
}

build_effdes <- function(nitems, mpoints, pplgrps, categos, groupvec){
  #builds treatment design structure for W
  #
  #mpoints>nitems>treat>catego
  #build group design
  tmp1 <- diag(pplgrps)
  tmp1[pplgrps, pplgrps] <- 0
  eff.tmp1 <- lapply(categos, function(x)(tmp1%x%x)) #list with categories per item, replicated per group
  eff.tmp2 <- as.matrix(bdiag(eff.tmp1))  #blockdiagonal with blocks equal to the categories
  eff.tmp3 <- diag(mpoints-1)%x%eff.tmp2  #blow up to mpoints
  nuller <- matrix(0, nrow=dim(eff.tmp2)[1], ncol=dim(eff.tmp3)[2]) #baseline (tp=1)
  gr.bu <- rbind(nuller, eff.tmp3) #combine baseline and effects
  #labelling of effects
  names1 <- unique(names(groupvec))
  #names1 <- paste0("G", pplgrps:1)
  names2 <- paste(names1, "I", sep=".")
  names3 <- paste0(names2, rep(1:nitems, each=pplgrps))
  names4 <- paste(names3, "t", sep=".")
  names5 <- paste0(names4, rep(2:mpoints, each=pplgrps*nitems))
  colnames(gr.bu) <- names5
  #columns with zeros (baseline group) are removed now
  rem.0 <- NA
  for(i in 1:dim(gr.bu)[2]) {rem.0[i] <- all(gr.bu[, i]==0)}
  gr.bu.red <- gr.bu[, which(rem.0==0)]
  return(gr.bu.red)
}


build_trdes <- function(nitems, mpoints, pplgrps, categos){
  #builds trend design structure for W
  #
  #mpoints>nitems>treat>catego
  tr.tmp1 <- lapply(categos, function(x) rep(x, pplgrps)) #replicate number of categories per item times the groups
  tr.tmp2 <- as.matrix(bdiag(tr.tmp1)) #build the blockdiaginal for all items
  tr.tmp3 <- diag(mpoints-1)%x%tr.tmp2 #blow it up to the time points necessary
  nuller <- matrix(0, nrow=dim(tr.tmp2)[1], ncol=dim(tr.tmp3)[2]) #baseline
  tr.bu <- rbind(nuller, tr.tmp3) #combine mpoints and baseline
  #structure: for each category multiply it with a vector of group indicators
  #hence the grouping is:
  #tau1 t2-t1, tau2 t2-t1, ..., tauk t2-t1, tau1 t3-t1, tau2 t3-t1, .. tauk t3-t1
 #cat("Design matrix columns are:","\n","tau_1^(t2-t1), tau_2^(t2-t1), ..., tau_k^(t2-t1), tau_1^(t3-t1), tau_2(t3-t1), ..., tau_k^(t3-t1), etc.","\n")
  #labeling
  names1 <- paste0("trend.I", 1:nitems)
  names2 <- paste(names1, "t", sep=".")
  names3 <- paste0(names2, rep(2:mpoints, each=nitems))
  colnames(tr.bu) <- names3
  return(tr.bu)
}

build_catdes <- function(nitems, mpoints, pplgrps, categos){
  #builds category design matrix
  #FIX ME: is a bit ugly, we might get the loops out somehow
  #
  #check if there are just binary items
  if(max(unlist(categos))<2) stop("items are (at most) binary and need no design")
  #currently equates cat.0 and cat.1
  warning("Currently c0 and c1 are equated for each item", "\n")
  max.all <- max(unlist(categos)) #maximum category number
  ls.ct.des <- list() #list of designs for each item
  #here we walk through each item and build up the category design
  for(i in 1:nitems){
    max.it <- sum(categos[[i]]!=0) #maximum category number of item i
    ct.des <- rbind(rep(0, dim(diag(max.all-1))[2]), diag(max.all-1)) #the design for the maximum number of categories in X
    rems <- max.all-max.it #the number of superfluous columns
    #here the superfluous columns are removed as the step from W to W*
    #the necessary rows with zeros however are maintained:
    #for a dichotomous item the structure is slightly different than for any other, since it returns an empty matrix of appropriate dimensions
    #for a polytomous item the superfluous columns are removed from the back
    ifelse(rems==max.all-1, ct.des<- as.matrix(ct.des[, -(1:max.all-1)]), ct.des<- as.matrix(ct.des[, 1:((max.all-1)-rems)]))
    ct.des.gr <- rep(1, pplgrps)%x%ct.des #blow it up to the number of groups
    ls.ct.des[[i]] <- ct.des.gr #list with all category designs for each item
  }
  ct.tmp2 <- as.matrix(bdiag(ls.ct.des)) #blockdiagonal matrix for a single mpoints
  ct.bu <- rep(1, mpoints)%x%ct.tmp2 #blow up to number of times points
  #try to first build first item, then second and so on, then blow up
  #labeling: pretty unelegant too
  names <- NA
  for(i in 1:nitems){
    cat <- max(categos[[i]])
    ifelse(cat==1, names1 <- "remove", names1 <- paste0("c", 2:cat))
    names2 <- paste0("I", i)
    names3 <- paste(names1, names2, sep=".")
    names<- c(names, names3)
  }
  names <- names[-1]
  if(length(grep("remove", names)>0)) names <- names[-grep("remove", names)]
  colnames(ct.bu) <- names
  return(ct.bu)
}
