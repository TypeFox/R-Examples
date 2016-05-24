

Hybrid.Tree=function(data,cluster0,model="nbinom"){

    n=data$Count
    s=data$Normalizer
    t=data$Treatment
    logFC=data$logFC
    n[n<=0]=1e-10
    nG=nrow(n)
    nb.disp=data$NB.Dispersion
    if( model=="poisson") nb.disp=rep(1e-10,nrow(n))
    nb.disp=matrix(rep(nb.disp,ncol(n)),nrow=nG)

   k=as.numeric(factor(cluster0))
   nK=length(unique(k))
   m=rep(0,nG)
   C=matrix(0,nK,length(unique(t)))
   rownames(C)=paste(1:nK)
   for(i in 1:nK){
      mc=cl.mb.est.mc(n[k==i,],s[k==i,],t,model,p=NULL,c=NULL,m=NULL,nb.disp=nb.disp[k==i,])
      m[k==i]=mc$m
      C[i,]=mc$c
   }
   tree=matrix(0,nK,3)
   colnames(tree)=c("k1","k2","distance")
   for(i in 1:(nK-1)){
     print(paste("level",i))
     D=dst.pairs(n,s,t,model,nb.disp,m,C,k,method="KL",nNearest=100)
     D=dst.pairs(n,s,t,model,nb.disp,m,C,k,method="Ward",nNearest=1,pairs=D[,1:2])
     tree[i,]=D
     k1=as.integer(D[1])
     k2=as.integer(D[2])
     k[k==k2]=k1
     mc=cl.mb.est.mc(n[k==k1,],s[k==k1,],t,model,p=NULL,c=NULL,m=m[k==k1],nb.disp=nb.disp[k==k1,])
     m[k==k1]=mc$m
     C[paste(k1),]=mc$c
     C=C[rownames(C)!=paste(k2),]
   }
  tree[nK,]=1
  return(tree)
}


#########################################################


sortNode=function(k1,k2){
 nK=length(k1)
 tree=as.list(1:nK)
 for(i in 1:(nK-1)){
   id1=in.tree(tree,k1[i])
   id2=in.tree(tree,k2[i])
   tree[[id1]]=c(tree[[id1]],tree[[id2]])
   tree=tree[-id2]
 }
 return(tree[[1]])
}
in.tree=function(tree,x){
  L=length(tree)
  id=0
  for(i in 1:L)
    if (x %in% tree[[i]]) id=i
  return(id)
}

move.two=function(s,s0,x,y){
  idx=s0==x
  idy=s0==y
  L=1:length(s)
  mdx=max(L[idx])
  s1=s[1:mdx]
  s2=s[idy]
  s3=s[L>mdx&(!idy)]
  s=c(s1,s2,s3)

  s1=s0[1:mdx]
  s2=s0[idy]
  s3=s0[L>mdx&(!idy)]
  s0=c(s1,s2,s3)
  s0[s0==y]=x
  return(list(s=s,s0=s0))
}

loc.node=function(k,k1,k2,h,num=NULL){
  nK=length(h)
  s=sortNode(k1,k2)
  k1X=k2X=1:nK
  x=cumsum(num[s])
  x=x/max(x)*nK
  xm=head(c(0,x),nK)/2+x/2
  for(i in 1:nK){
   k1X[i]=xm[s==k1[i]]
   k2X[i]=xm[s==k2[i]]
  }
  k1H=k2H=rep(0,nK)
  H=matrix(0,nK,nK)
for( i in 2:nK){
 id=1:(i-1)
 id=tail(id[k1[1:(i-1)]==k1[i]],1)
 if(length(id)>0){
   k1H[i]=k1H[id]+k2H[id]+h[id]
   k1X[i]=k1X[id]/2+k2X[id]/2
 }
 id=1:(i-1)
 id=tail(id[k1[1:(i-1)]==k2[i]],1)
 if(length(id)>0){
   k2H[i]=k1H[id]+k2H[id]+h[id]
 #  k2H[i]=k1H[id]+k2H[id]
   k2X[i]=k1X[id]/2+k2X[id]/2
  }
 }
 return(list(k1h=k1H,k2h=k2H,k1x=k1X,k2x=k2X))
}





plotbr=function(cl,cr,ch){
  sx0=c(cl[1],cr[1],cl[1])
  sy0=c(cl[2],cr[2],ch)
  sx1=c(cl[1],cr[1],cr[1])
  sy1=c(ch,ch,ch)
  segments(sx0,sy0,sx1,sy1)
}

##############################################
###### pt: pattern; h: height; merge b to a
##############################################
leaf.color=function(centers){
  c=c(centers)
  c=rank(c)
  c=c/max(c)
  c=1-.9*c
  c=matrix(c,nrow=nrow(centers),ncol=ncol(centers))
  return(c)
}


plotHybrid.Tree=function(merge,cluster,logFC,tree.title=NULL,colorful=FALSE){
 k1=merge[,1]
 k2=merge[,2]
 if(is.null(tree.title)) tree.title=""
 h=rep(1e-3,length(k1))
 num=table(cluster)
 centers=logFC
 nT=ncol(centers)
 nJ=nrow(centers)
 nK=K=length(k1)
 k=1:nK
 s=sortNode(k1,k2)
 foo=loc.node(k,k1,k2,h,num)
 k1h=foo$k1h
 k2h=foo$k2h
 k1x=foo$k1x
 k2x=foo$k2x
 
 maxH=max(c(k1h+k2h)[1:(K-1)])
 pt=leaf.color(centers)
 PT=c()
 for(i in 1:nK) PT=rbind(PT,pt[cluster==s[i],])
 pt=PT
 par(mar=c(1,1,1,1))
 plot(0,0,xlim=c(0,nK),ylim=c(-1.25*maxH/2,1.1*maxH),cex=0,xaxt="n",yaxt="n",xlab="",ylab="",main="",cex.main=1)
 text(nK/2,1.1*maxH,tree.title,cex=1.2)
 text(-6,-((1:nT)-.4)*maxH/nT/2,1:nT,pos=4,cex=1)
 for( i in 1:(K-1)){
   clx=k1x[i]
   crx=k2x[i]
   cly=k1h[i]
   cry=k2h[i]
   ch=h[i]+k1h[i]+k2h[i]
   plotbr(c(clx,cly),c(crx,cry),ch)
 }

 x=1:nrow(pt)
 x=x/max(x)*K
if(!colorful) image(c(0,x),seq(-nT,0,1)*maxH/nT/2-.02*maxH,pt,col=grey((0:nrow(pt))/nrow(pt)),add=TRUE)
if(colorful) image(c(0,x),seq(-nT,0,1)*maxH/nT/2-.02*maxH,pt,col=heat.colors(nrow(pt)),add=TRUE)

 x=cumsum(num[s])
 x=x/max(x)*K
 xm=head(c(0,x),K)/2+x/2
 ym=(-1.22+.1*rep(c(-1,0,1,0),K)[1:nK])*maxH/2
 text(xm,ym,s,cex=.5,srt=90,pos=2,offset=0)
 segments(c(0,x),-1.03*maxH/2,c(0,x),-1.08*maxH/2)

if(!colorful) image(c(0,nK),c(-.03*maxH/2,0),matrix(.5,1,1),col=grey(c(0:16)/16),add=TRUE)
if(colorful) image(c(0,nK),c(-.03*maxH/2,0),matrix(.5,1,1),col=heat.colors(16),add=TRUE)
 segments(c(0,x),-.03*maxH/2,c(0,x),0)
 segments(0,0,nK,0)
 ######## legend
 ct=c(centers)
 cl=leaf.color(centers)
 ord=order(ct)
 ct=ct[ord]
 cl=cl[ord]
 nbar=12
 
 id=round(seq(1,length(ct),length.out=nbar+2))
 id=id[-1]
 id=id[1:nbar]
 ct=ct[id]
 cl=cl[id]
 x=c(.1*K,.13*K)
 y=seq(.5*maxH,1.1*maxH,length.out=nbar+1)
if(!colorful) image(x,y,matrix(cl,1,nbar),col=grey((0:32)/32),add=TRUE )
if(colorful) image(x,y,matrix(cl,1,nbar),col=heat.colors(32),add=TRUE )

     tx=paste(round(ct,2)+1e-4*sign(round(ct,2)))
     tx=unlist(strsplit(tx,"\\."))
     tx1=tx[(1:nbar)*2-1]
     tx2=tx[(1:nbar)*2]
     tx2=substr(tx2,1,2)
     tx=paste(tx1,tx2,sep=".")
     tx[1]=paste("<",tx[1])
     tx[nbar]=paste(">",tx[nbar])
 text(.09*nK,y[1:nbar]/2+y[-1]/2,tx,pos=2,cex=.6)
}

