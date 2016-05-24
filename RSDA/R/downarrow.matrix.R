downarrow.matrix <-
function(sym.data) {
  k<-max(sym.data$sym.var.length)
  nn<-sym.data$N
  mm<-sym.data$M
  n.sym.objects<-k*nn
  n.sym.var<-mm
  meta.M<-3*mm
  sym.var.types<-rep("$I",n.sym.var)
  sym.var.length<-rep(2,n.sym.var)
  sym.var.names<-sym.data$sym.var.names
  sym.var.starts<-seq(2,meta.M,3)
  sym.obj.names<-seq(1,n.sym.objects)  
  m1<-as.data.frame(rep(0,n.sym.objects))
  for(j in 1:mm) {
    if(sym.var(sym.data,j)$var.type=='$C') {
      pos<-1
      c1<-rep("$I",n.sym.objects)
      c2<-rep(0,n.sym.objects)
      c3<-rep(0,n.sym.objects)
      for(i in 1:nn) {
        for(s in 1:k) {
          c2[pos]<-sym.var(sym.data,j)$var.data.vector[i]
          c3[pos]<-sym.var(sym.data,j)$var.data.vector[i]
          pos<-pos+1
        }
      }
      m1<-cbind(m1,c1,as.numeric(c2),as.numeric(c3))
    }  
    if(sym.var(sym.data,j)$var.type=='$I') {
      pos<-1
      c1<-rep("$I",n.sym.objects)
      c2<-rep(0,n.sym.objects)
      c3<-rep(0,n.sym.objects)
      for(i in 1:nn) {
        for(s in 1:k) {
          c2[pos]<-sym.var(sym.data,j)$var.data.vector[i,1]
          c3[pos]<-sym.var(sym.data,j)$var.data.vector[i,2]
          pos<-pos+1
        }
      }
      m1<-cbind(m1,c1,as.numeric(c2),as.numeric(c3))
    }
    if(sym.var(sym.data,j)$var.type=='$H') {
      pos<-1
      c1<-rep("$I",n.sym.objects)
      c2<-rep(0,n.sym.objects)
      c3<-rep(0,n.sym.objects)
      ff<-sym.data$sym.var.length[j]
      for(i in 1:nn) {
        prob.acum<-0
        for(s in 1:k) {
          if(s<=ff) {
            prob.acum<-prob.acum+sym.var(sym.data,j)$var.data.vector[i,s]
            c3[pos]<-prob.acum
          }
          else {
            c3[pos]<-prob.acum
          }  
          pos<-pos+1  
        }
      }
      m1<-cbind(m1,c1,as.numeric(c2),as.numeric(c3))
    }
  }  
  meta.data<-m1[,-1]
  vnames<-rep("$I",meta.M)  
  pos<-2
  for(i in 1:mm) {
    vnames[pos]<-sym.data$sym.var.names[i]
    vnames[pos+1]<-sym.data$sym.var.names[i]
    pos<-pos+3
  }
  colnames(meta.data)<-vnames 
  data.matrix<-meta.data
  pos<-1
  for(i in 1:mm) {
     data.matrix<-data.matrix[,-pos]
     pos<-pos+2
  }  
  sym.data<-list(N=n.sym.objects,M=n.sym.var,sym.obj.names=sym.obj.names,
                 sym.var.names=sym.var.names,sym.var.types=sym.var.types,
                 sym.var.length=sym.var.length,sym.var.starts=sym.var.starts,
                 meta=meta.data,data=data.matrix)
  return(sym.data)  
}
