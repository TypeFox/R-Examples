bootdegK <- function(sam.out, num.sam, n.boot) {
  n.seeds<-sam.out$n.seeds
  n.neigh<-sam.out$n.neigh
  if(length(num.sam)==1) num.sam<-1:num.sam

  empd<-as.list(rep(NA,length(num.sam)))
  i<-1
  for(m in num.sam){


    val.seeds<-sam.out$val.seeds[[m]]
    val.nonseed<-sam.out$val.nonseed[[m]]
    freq.deg.seeds<-sam.out$samples[[m]]$freq.deg.seeds
    freq.deg.nonseed<-sam.out$samples[[m]]$freq.deg.nonseed

    bsam.seeds<-myBsample(val.seeds, n.seeds,n.boot, prob=freq.deg.seeds)
    # ^matrix n.boot x n.seeds
    bsam.nonseed.nw<-myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob=freq.deg.nonseed)
    # ^matrix n.boot x sum(freq.deg.nonseed)
    bsam.nonseed.w<-myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob=freq.deg.nonseed/val.nonseed)
    # ^matrix

    p0.B<-rep(0,n.boot)
    if(any(val.seeds==0)){
      # ^if any seeds has degree zero
      p0.B<-rowSums(bsam.seeds==0)/n.seeds
      # ^the estimation from the bootstrap samples
    }
    p0.real<-sam.out$p0.real
    p0.seeds<-sam.out$p0.seeds[[m]]

    values<-sam.out$values[[m]]
    # ^all the possible degree values to resample

    Fseed<-t(apply(bsam.seeds,1,table.row,vect=values))
    # ^frequency (sorted according to values)

    if (is.null(bsam.nonseed.nw)){
      # browser()
      Fnonseed.nw<-0
    }else{Fnonseed.nw<-t(apply(as.matrix(bsam.nonseed.nw),1,table.row,vect=values))
    }
    if (is.null(bsam.nonseed.w)){
      Fnonseed.n<-0
    }else{Fnonseed.w<-t(apply(as.matrix(bsam.nonseed.w),1,table.row,vect=values))
    }

    ####################################################################
    #####       combining information from seeds and nonseeds     ######
    ####################################################################

    #mean degree computed from the original sampled seeds:
    ekseed<-sam.out$ekseed[[m]]

    colzero<-NULL
    if(any(values==0)){
      colzero<-which(values==0)
      vals<-values[-colzero]
      f.seeds<-Fseed[,-colzero]
      f.nonseed.nw<-Fnonseed.nw[,-colzero]
      f.nonseed.w<-Fnonseed.w[,-colzero]
    }else{
      vals<-values
      f.seeds<-Fseed
      f.nonseed.nw<-Fnonseed.nw
      f.nonseed.w<-Fnonseed.w
    }
    #browser()
    ################################################
    # WB  #   seeds and weighted nonseeds     ######
    ################################################
    ### consider the p0 fixed from the seeds information:
    #empd.w.p0s<-(f.seeds+(1-p0.seeds)*f.nonseed.w)/(n.seeds+sum(freq.deg.nonseed))
    empd.w.p0s<-(f.seeds+f.nonseed.w*(1-p0.B))/(n.seeds+sum(freq.deg.nonseed))
    #### consider the p0 fixed from the real information:
    #empd.w.p0r<-(f.seeds+(1-p0.real)*f.nonseed.w)/(n.seeds+sum(freq.deg.nonseed))
    ################################################
    # NWB #  seeds and non weighted nonseeds    ####
    ################################################
    #########################################
    #p0 estimated from orginal sampled seeds#
    #########################################
    #E(K) estimated from bootstrap samples from the seeds
    #empd.nw.p0sEkb<-(f.seeds+(1-p0.seeds)*apply(bsam.seeds,1,FUN=mean)*t(t(f.nonseed.nw)/vals))/(n.seeds+ apply(bsam.seeds,1,FUN=mean)*rowSums(t(t(f.nonseed.nw)/vals)))
    empd.nw.p0sEkb<-(f.seeds+t(t(f.nonseed.nw)/vals)*(1-p0.B)*apply(bsam.seeds,1,FUN=mean))/(n.seeds+rowSums(t(t(f.nonseed.nw)/vals))*apply(bsam.seeds,1,FUN=mean))
    #E(K) estimated from the original seeds sample
    #empd.nw.p0sEks<-(f.seeds+(1-p0.seeds)*ekseed*t(t(f.nonseed.nw)/vals))/(n.seeds+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))
    empd.nw.p0sEks<-(f.seeds+ekseed*t(t(f.nonseed.nw)/vals)*(1-p0.B))/(n.seeds+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))
    #########################################
    #p0 taken as known                      #
    #########################################
    #E(K) estimated from bootstrap samples from the seeds
    #empd.nw.p0rEkb<-(f.seeds+(1-p0.real)*apply(bsam.seeds,1,FUN=mean)*t(t(f.nonseed.nw)/vals))/ (n.seeds+apply(bsam.seeds,1,FUN=mean)*rowSums(t(t(f.nonseed.nw)/vals)))
    #E(K) estimated from the original seeds sample
    #empd.nw.p0rEks<-(f.seeds+(1-p0.real)*ekseed*t(t(f.nonseed.nw)/vals))/(n.seeds+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))

    if(any(values==0)){
      empd.w.p0s<-cbind("0"=p0.B,empd.w.p0s)
      empd.nw.p0sEkb<-cbind("0"=p0.B,empd.nw.p0sEkb)
      empd.nw.p0sEks<-cbind("0"=p0.B,empd.nw.p0sEks)
    }
    empd[[i]]<-list(empd.w.p0s=empd.w.p0s,empd.nw.p0sEkb=empd.nw.p0sEkb, empd.nw.p0sEks= empd.nw.p0sEks)
    i<-i+1
  }  # for(m in num.sam)
  list(values = sam.out$values[num.sam], empd = empd, num.sam = num.sam, n.boot = n.boot,
       n.neigh = n.neigh, seeds1 = sam.out$seeds1[num.sam,],
       nodes_of_LSMI = sam.out$nodes_of_LSMI[num.sam])
}
