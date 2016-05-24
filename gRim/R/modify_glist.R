##################################################################
####
#### Update generator list by adding/deleting edges and terms
####
#### FIXME: Perhaps add... should check if ee/term is in the list
#### in which case a special value should be returned
####
##################################################################

### Updates an entire glist with the elements (edges, terms) in
### items
###
modify_glist <- function(glist, items, details=0){
  glist   <- lapply(glist, as.character)
  items   <- .parse.change.list(items, details)
  action  <- names(items)
  for (ii in seq_along(items)){
    curr.action  <- action[ii]
    curr.item    <- items[[ii]]
    glist        <- .modify_glistPrim(glist, curr.action, curr.item, details)
  }
  glist
}

### Updates a glist (generating class) with the elements in curr.item. These can be of
### the type curr.action where valid choices are add.edge, drop.edge, add.term and drop.term
###
.modify_glistPrim <- function(glist, curr.action, curr.item, details=0){
  fname <- paste(".",curr.action,"_glist",sep="")
  #cat(sprintf("fname: %s\n", fname))
  .infoPrint(details,1,cat(sprintf("action: %s \n", curr.action)))

  for (kk in seq_along(curr.item)){
    curr <- curr.item[[kk]]
    ##cat(sprintf("action: %s item: %s\n", fname, paste(curr, collapse=" ")))
    glist <- do.call(fname, list(glist, curr))
  }
  return(glist)
}

.add.edge_glist <- function(glist, ee){
  extra <- list()
  cnt <- 1
  for (ii in seq_along(glist)){
    if (ee[1] %in% glist[[ii]]){
      for (jj in seq_along(glist)){
        if (ee[2] %in% glist[[jj]]){
           extra[[cnt]] <- uniquePrim(c(ee, intersectPrim(glist[[ii]], glist[[jj]])))
           cnt <- cnt + 1
        }
      }
    }
  }

  glist.new <- removeRedundant(c(glist, extra))
  return(glist.new)
}

.drop.edge_glist <- function(glist, ee){
  .drop.term_glist(glist, ee)
}

.add.term_glist <- function(glist, term){
  if (isin(glist,term))
    return(glist)
  else
    return(removeRedundant(c(list(term), glist)))
}

.drop.term_glist <- function(glist, term){
    extra <- list()
    cnt <- 1
    changed <- rep(0, length(glist))

    for (ii in seq_along(glist)){
        gg <- glist[[ii]]
        ##cat("term:\n"); print(term);
        ##cat("gg:\n"); print(gg)
        if (subsetof(term,gg)){
            ##cat("term is subset of gg...\n")
            changed[ii] <- 1
            short1 <- combnPrim(gg, length(gg)-1, simplify=FALSE)
            if (length(term)==length(gg)){
                extra[[cnt]] <- short1
            } else {
                keep   <- unlistPrim(lapply(short1, function(sss) !subsetof(term,sss)))
                                        #print(keep)
                short1 <- short1[keep]
                extra[[cnt]] <- short1
            }
            cnt <- cnt + 1
        }
    }

    glist.new <- c(glist[changed==0], unlist(extra, recursive=FALSE,use.names=FALSE))
    glist.new <- removeRedundant(glist.new)
    return(glist.new)
}

.aedge_glist <- .add.edge_glist
.dedge_glist <- .drop.edge_glist
.aterm_glist <- .add.term_glist
.dterm_glist <- .drop.term_glist


### An ad.list can have elements with names add.edge, drop.edge, add.term and drop.term
### These can be formulae, and .parse.change.list will transform these into appropritate
### lists.
###
.parse.change.list <- function(items,details=0){
                                        #cat("In function: .parse.change.list:\n")
  .foo <- function(curr.action, curr.item){
    switch(curr.action,
           "add.edge"=,"drop.edge"=,"aedge"=,"dedge"={
             zzz <- unlist(lapply(
                                  rhsf2list(curr.item),
                                  names2pairs),recursive=FALSE)
           },
           "add.term"=,"drop.term"=,"aterm"=,"dterm"={
             zzz <- rhsf2list(curr.item)
           })
    zzz
  }
  nam <- names(items)
  valid <- c("add.edge","drop.edge","add.term","drop.term",
             "aedge","dedge","aterm","dterm")

  for (ii in 1:length(items)){
    curr.action <- nam[ii]
    aaa <- matchPrim(curr.action, valid)
    if (is.na(aaa))
      stop(sprintf("Item %i has name '%s' which is not valid\n",ii, curr.action))
    curr.item <- items[[ii]]
    .infoPrint(details,1, cat(sprintf("parsing action %s on item %s\n", curr.action, toString(curr.item))))
    items[[ii]] <- .foo(curr.action, curr.item)
  }
                                        #cat("On exit:\n"); print(items)
  return(items)
}









