.onUnload <- function (libpath) {
    library.dynam.unload("rankdist", libpath)
}

paramTow = function(param.true){
    w.true = rev(cumsum(rev(param.true)))
    w.true
}

wToparam = function(w.true){
    param.true = numeric(length(w.true))
    param.true[1:(length(w.true)-1)] = -diff(w.true)
    param.true[length(param.true)] = w.true[length(w.true)]
    param.true
}

UpdateCount = function(dat,count){
  dat@count = count
  dat@nobs = sum(count)
  if (length(dat@topq)!=1 && min(dat@topq) < dat@nobj-1){
    dat@subobs = numeric(length(dat@topq))
    for (i in 1:length(dat@topq)){
      dat@subobs[i] = sum(dat@count[ dat@q_ind[i]: (dat@q_ind[i+1]-1) ])
    }
  }
  dat
}



# used in SearchPi0: make a optimization result into a model
AddInfo=function(solveres,dat,pi0){
    solveres$nobs = dat@nobs
    solveres$nobj = dat@nobj
    solveres$pi0.ranking = pi0
    solveres
}

# neighbour for incomplete rankings
SearchPi0=function(dat,init,ctrl){
#     if (class(ctrl)=="RankControlWtau"){
#       mod <- SingleClusterModel(dat,init,ctrl,0)
#       return(mod)
#     }
    
    n = dat@nobj
    curr_best_ranking = init@modal_ranking.init[[1]] 
    if (ctrl@SearchPi0_show_message){
      message("<<< initial ranking ",curr_best_ranking," >>>")
    }
    if (max(dat@topq) < n-1){
        curr_best_ranking[curr_best_ranking>max(dat@topq)+1]=max(dat@topq)+1
    }
    curr_solve <- SingleClusterModel(dat,init,ctrl,curr_best_ranking)
    curr_model = AddInfo(curr_solve,dat,curr_best_ranking)
    FUN = ctrl@SearchPi0_FUN
    curr_goodness = FUN(curr_model)
    hashtable = hash::hash()
    hash::.set(hashtable,keys = RanktoHash(curr_best_ranking),values=TRUE)
    SearchPi0_step = 0
    while(TRUE){
        SearchPi0_step = SearchPi0_step+1
        if (SearchPi0_step > ctrl@SearchPi0_limit){
            if (ctrl@SearchPi0_show_message){
                message("Search Pi0 limit has been reached. Stop at current best: ",this_ranking)
            }
            break
        }
        if (ctrl@SearchPi0_neighbour=="Cayley"){
            neighbours = CayleyNeighbour(curr_best_ranking)
        } else {
            neighbours = KendallNeighbour(curr_best_ranking)
        }
        testkeys = RanktoHash(neighbours)
        tested = hash::has.key(testkeys,hashtable)
        if (all(tested)) break
        hash::.set(hashtable,keys=testkeys[!tested],values=rep(TRUE,length(testkeys[!tested])))
        for (i in 1:nrow(neighbours)){
            # tested neighbours cannot be better
            if (tested[i]) next
            this_ranking = neighbours[i,]
            if (ctrl@SearchPi0_show_message){
                message("\tNow Checking Neighbour ",this_ranking)
            }
            this_solve <- SingleClusterModel(dat,init,ctrl,this_ranking)
            this_model = AddInfo(this_solve,dat,this_ranking)
            this_goodness = FUN(this_model)
            if (this_goodness > curr_goodness){
                curr_goodness = this_goodness
                curr_best_ranking = this_ranking
                curr_model = this_model
                if (ctrl@SearchPi0_show_message){
                    message("***Best changed to ",curr_best_ranking,"***")
                }
                if (ctrl@SearchPi0_fast_traversal)
                    break
            }
        }
        if (ctrl@SearchPi0_show_message){
            message("--> Moved to ",curr_best_ranking," <--")
        }
    }
    curr_model$SearchPi0_step = SearchPi0_step 
    curr_model
}

# TODO need to change: does not work for d==1
t.gen = function(d){
    t.lst = list()
    t.lst[[d]] = matrix(rep(1:d,d),ncol = d, nrow = d,byrow=T)
    left.mask = matrix(rep(0,d^2),ncol = d, nrow = d)
    left.mask[2:d,1:(d-1)] = diag(rep(1,d-1))
    t.lst[[d]][upper.tri(left.mask)] = 0
    for ( i in 1:(d-1)){
        t.lst[[d-i]] = left.mask%*%t.lst[[d-i+1]]
        diag(t.lst[[d-i]]) = c(rep(0,i),1:(d-i))
        t.lst[[d-i]][upper.tri(left.mask)] = 0
    }
    t.lst
}
GHC = function(param,t.lst){
    d = length(param) # d = t - 1
    K = matrix(rep(0,d^2),ncol = d, nrow = d)
    for ( i in 1:d){
        K = -1 * param[i] * t.lst[[i]] + K
    }
    K = exp(K)
    K[upper.tri(K)] = 0
    gradiant = numeric(d)
    ones = rep(1,d)
    denom = rowSums(K) + ones
    B = matrix(ncol=d,nrow=d)
    for (i in 1:d){
        B[,i] = rowSums(-1 * K * t.lst[[i]])
    }
    for ( i in 1:d){
        gradiant[i] = sum(B[,i] / denom)
    }
    gradiant
}



# find the weighted kendall distance between p1 and p2
# p1 and p2 are orderings
KwDist = function(p1, p2,w){
    n = length(p1)
    distance = 0
    for (i in p2){
        pos1 = which(p1 == i)
        pos2 = which(p2 == i)
        relative_pos1 = (1:n - pos1)[order(p1)]
        relative_pos2 = (1:n - pos2)[order(p2)]
        Ji = which(relative_pos1 * relative_pos2 < 0)
        Ii = length(Ji)
        Li = (pos1 + pos2 + Ii)/2
        c1 = ifelse(pos1<=(Li-1), sum(w[pos1:(Li-1)]),0)
        c2 = ifelse(pos2<=(Li-1), sum(w[pos2:(Li-1)]),0)
        distance = distance + (c1 + c2)/2
    }
    distance
}

BreakTieEqualProb <- function(dat){
    ind_comp <- which(dat@topq == dat@nobj-1)
    if (max(dat@topq) == dat@nobj-1){
        ind_comp_start <- dat@q_ind[ind_comp]
        ind_comp_end <- dat@q_ind[ind_comp + 1] - 1
        comp_ranking <- dat@ranking[ind_comp_start:ind_comp_end, ]
        comp_count <- dat@count[ind_comp_start:ind_comp_end]
    } else {
        comp_ranking <- permute::allPerms(dat@nobj)
        comp_ranking <- rbind(1:dat@nobj, comp_ranking)
        comp_count <- rep(0, nrow(comp_ranking))
    }
    comp_hash <- RanktoHash(comp_ranking)
    
    for (i in 1:length(dat@topq)){
        if(dat@topq[i] == dat@nobj-1) 
            next
        ind_start <- dat@q_ind[i]
        ind_end <- dat@q_ind[i+1] - 1
        this_q <- dat@topq[i]
        this_inc <- 1/factorial(dat@nobj - this_q)
        # generate permutations for tied group
        tie_perm <- permute::allPerms((this_q+1):dat@nobj) + this_q
        tie_perm <- rbind(tie_perm, (this_q+1):dat@nobj)
        # iterate through top-q rankings
        for (this_partial_ind in ind_start:ind_end){
            this_partial <- dat@ranking[this_partial_ind, ]
            this_count <- dat@count[this_partial_ind]
            ind_tie <- which(this_partial == this_q + 1)
            # iterate through possible tie-breakings
            for (ind_break in 1:nrow(tie_perm)){
                copy_partial <- this_partial
                this_break <- tie_perm[ind_break, ]
                ptr_break <- 1
                # iterate through tied positions
                for (this_tie_ind in ind_tie){
                    copy_partial[this_tie_ind] = this_break[ptr_break]
                    ptr_break <- ptr_break + 1
                }
                this_hash <- rankdist::RanktoHash(copy_partial)
                ind_incre <- which(comp_hash == this_hash)
                comp_count[ind_incre] = comp_count[ind_incre] + this_inc*this_count
            }
        }
    }
    # handle complete rankings
    ind_nonempty_count = which(comp_count != 0)
    comp_count = comp_count[comp_count > 0]
    comp_ranking = comp_ranking[ind_nonempty_count, ]
    comp_dat <- new("RankData", ranking=comp_ranking, count=comp_count)
}

