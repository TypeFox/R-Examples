build_W <- function(X,nitems,mpoints,grp_n,groupvec,itmgrps)
  {
     if(missing(grp_n)) grp_n<- table(groupvec) 
     if(!is.numeric(grp_n)) stop("Please specify the number of subjects per group.")
     if(missing(nitems))stop("Please specify the number of items.")
     if(any(grp_n==0)) stop("There are groups with zero sample size.")
     if(missing(mpoints)) stop("Please specify the number of time points. If there are none, you might want to use PCM() or LPCM().")
    pplgrps <- length(grp_n)
    #builds the LLRA design matrix from scratch
    categos <- get_item_cats(X,nitems,grp_n)
    #trend effects design
    tr.des <- build_trdes(nitems,mpoints,pplgrps,categos)
    #tretment effects design
    gr.des <- build_effdes(nitems,mpoints,pplgrps,categos,groupvec)
    #category design
    if(length(unique(unlist(categos)))==1&&sum(unique(unlist(categos)))==1) return(cbind(gr.des,tr.des))
    ct.des <- build_catdes(nitems,mpoints,pplgrps,categos)
    #all together now!
    des <-cbind(gr.des,tr.des,ct.des)
    des
  }
