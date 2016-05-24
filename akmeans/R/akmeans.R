akmeans <-
function(x, ths1=0.2, ths2=0.2, ths3=0.7, ths4=0.2, min.k=5, max.k=100, iter.max=100, nstart=1, mode=1, d.metric=1, verbose=TRUE){
## x: data matrix n by p: all elements should be numeric
## ths1: threshold to decide whether to increase k or not: check sum((sample-assigned center)^2) < ths1*sum(assigned center^2)
## ths2: threshold to decide whether to increase k or not: check all components of |sample-assigned center| < ths2
## ths3: threshold to decide whether to increase k or not: check inner product of (sample,assigned center) > ths3 , this is only for cosine distance metric 
## ths4: threshold to decide whether to increase k or not: check sum(abs(sample-assigned center)) < ths4 
## mode: 1: use ths1, 2: use ths2, 3: use ths3, 4: use ths4
## d.metric: 1: use euclidean distance metric, otherwise use cosine distance metric
## iter.max, nstart: will be delivered to kmeans
## min.k: minimum k=> start kmeans from min.k
## max.k: maximum k=> will stop at max.k
  
## algorithm is quite simple
## 1. Set min.k and max.k.
## 2. Run K-means with K = min.k
## 3. For each cluster, check the threshold condition.
## 4. If all clusters satisfy the threshold condition => Done, return the result
## 5. Check K>max.k => If yes, stop. If no, go to step 5.
## 6. For any cluster violating the threshold condition, run K'-means with K'=2 on those cluster members, 
## which means K will increase by the number of violating clusters.
## 7. Run K-means setting the present cluster centers as the initial centers and go to step 4.

  require('stats')

  if (!(mode%in%c(1:4))) stop('Choose a valid mode')
  now.k = min.k; n = nrow(x); p = ncol(x)

  if (d.metric==1){
    res = kmeans(x,centers=now.k,iter.max,nstart)
  } else {
    x = quick.norm(x,mod=1)
    res = norm.sim.ksc(x,now.k,iter.max=iter.max)
  }
  now.centers = res$centers

  while(1){ ## loop till converge or max.k
    if (mode==1){ ## use ths1
      vio.idx = which(apply((x-res$centers[res$cluster,])^2,1,sum) > ths1*apply(res$centers[res$cluster,]^2,1,sum))
    } else if (mode==2){ ## use ths2
      vio.idx = which(apply(abs(x-res$centers[res$cluster,]),1,function(i){any(i>ths2)})==TRUE)
    } else if (mode==4){ ## use ths4
      vio.idx = which(apply(abs(x-res$centers[res$cluster,]),1,sum)>ths4)
    } else { ## use ths3 this is especially for cosine distance
      vio.idx = which(apply(x*res$centers[res$cluster,],1,sum) < ths3)
    }

    if (length(vio.idx)>0){ ## some cluster doesn't satisfy the threshold condition
      re.cluster.list = unique(res$cluster[vio.idx]) ## cluster list to divide more
      if (verbose) print(paste('# of clusters violating given threshold condition: ',length(re.cluster.list)))
      if (now.k+length(re.cluster.list)>max.k) stop(paste('not converged till max.k. it needs at least ',length(re.cluster.list),' more clusters'))

      now.k = now.k+length(re.cluster.list)
      add.centers = matrix(0,2*length(re.cluster.list),p)
      if (verbose) print(paste('now.k=',now.k))
      cnt=0
      for (i in re.cluster.list){
        cnt = cnt + 1
        if (d.metric==1){
          if (sum(res$cluster==i)>2){
            add.centers[(2*cnt-1):(2*cnt),] = kmeans(x[res$cluster==i,],2,iter.max,nstart)$centers
          } else { ## kmeans can't run 2-means on 2 samples
            add.centers[(2*cnt-1):(2*cnt),] = x[res$cluster==i,]
          }
        } else {
          add.centers[(2*cnt-1):(2*cnt),] = norm.sim.ksc(x[res$cluster==i,],2,iter.max=iter.max)$centers
        } 
      }
      ## with updated centers, run kmeans again
      if (d.metric==1){
        res = kmeans(x,centers=rbind(now.centers[-re.cluster.list,],add.centers),iter.max)
      } else {
        res = norm.sim.ksc(x,now.k, init.cen=rbind(now.centers[-re.cluster.list,],add.centers),iter.max = iter.max)
      }
      now.centers = res$centers
    } else { ## converged
      if (verbose) print(paste('converged at k=',now.k))
      return(res)
    }  
  }
}

