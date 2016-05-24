global.efficiency <- function(adj.mat,weight.mat){
  
  n.regions<-dim(adj.mat)[1]
  if(sum(adj.mat)==0){
    glob.eff<-0
    ideal.eff<-1
    lp<-rep(0,n.regions)
    new.lp<-rep(0,n.regions)
  }
  else
  {
    
    if(dim(weight.mat)[1]!=n.regions) stop("Problem, the weight matrix and adjacency matrix must have the same size")
    if(dim(weight.mat)[2]!=n.regions) stop("Problem, the weight matrix and adjacency matrix must have the same size")
    z<-.C("Refficiencyfun",
        as.integer(n.regions),
        as.integer(adj.mat),
        as.double(weight.mat),
        Lp=double(n.regions),
        eff=double(1),
        matlp=double(n.regions*n.regions), PACKAGE="brainwaver")
    lp<-z$Lp
    matlp<-z$matlp

    new.lp<-rep(0,n.regions)
    for(fun in 1:n.regions){
      d<-matlp[(n.regions*(fun-1)+1):(n.regions*(fun))]
      d<-d[-fun]
      if(sum(d!=1e5)==0){
        new.lp[fun]<-0
      }else{
        new.lp[fun]<-(sum(1/d[d!=1e5])/(n.regions-1))
      }
    }

    glob.eff<-z$eff
    ideal.eff<-0
    for(i in 1:n.regions){
      for(j in 1:n.regions){
        if(i!=j){
          tmp<-weight.mat[i,j]
          add<-1/tmp
          ideal.eff<-ideal.eff+add
        }
      }
    }

    ideal.eff<-ideal.eff/(n.regions*(n.regions-1))
  }
  list(nodal.eff=new.lp,eff=(glob.eff/ideal.eff),new.lp=new.lp)
}


local.efficiency<- function(adj.mat,weight.mat){
  n.regions<-dim(adj.mat)[1]
  loc.eff<-rep(0,n.regions)
  if(sum(adj.mat)==0){
    eff<-0
  }else{
    eff<-0
    for(i in 1:n.regions){
      n.nodes<-sum(adj.mat[i,])
      if((n.nodes==0)||(n.nodes==1)){
        tmp.eff<-0
      }else{
        num.nodes<-rep(0,n.nodes)
        count<-1
        for(j in 1:n.regions){
          if(adj.mat[i,j]==1){
            num.nodes[count]<-j
            count<-count+1
          }
        }
        subgraph<-matrix(0,n.nodes,n.nodes)
        subweight<-matrix(0,n.nodes,n.nodes)
        count<-1
        for(j in 1:(n.nodes-1)){
          tmp<-num.nodes[j]
          count<-j+1
          for(k in (j+1):n.nodes){
            tmp1<-num.nodes[count]
            if(adj.mat[tmp,tmp1]==1){
              subgraph[j,k]<-1
            }
            subweight[j,k]<-weight.mat[tmp,tmp1]
            count<-count+1
          }
        }
        subgraph<-subgraph+t(subgraph)
        subweight<-subweight+t(subweight)
        tmp.eff.bis<-global.efficiency(subgraph,subweight)
        tmp.eff<-tmp.eff.bis$eff
        loc.eff[i]<-tmp.eff
      }
      eff<-eff+tmp.eff
    }
  }
  list(eff=(eff/n.regions),loc.eff=loc.eff)
}

cost.evaluator  <- function(x){
  return(x)
}

global.cost<-function(adj.mat,weight.mat){
  if(sum(adj.mat)==0){
    tmp<-0
    tmp1<-1
  }else{
    n.regions<-dim(adj.mat)[1]
    tmp<-0
    tmp1<-0
    for(i in 1:n.regions){
      for(j in 1:n.regions){
        if(i!=j){
          tmp<-tmp+adj.mat[i,j]*cost.evaluator(weight.mat[i,j])
          tmp1<-tmp1+cost.evaluator(weight.mat[i,j])
        }
      }
    }
  }
  return(tmp/tmp1)
}
