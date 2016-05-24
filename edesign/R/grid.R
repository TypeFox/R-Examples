sample.grid <- function(x,y,stepx=NULL,stepy=NULL,method="q",gp=NULL)
  {

    x <- as.double(x)
    y <- as.double(y)
    gp <- as.integer(gp)
    if(length(x)!=length(y))
      {
        stop("x and y have to be the same length")
      }
    testmethod <- switch(method, q="ok",r="ok",t="ok",h="ok",c="ok","error")
    if(testmethod=="error")
      stop("wrong argument for \"method\", should be one of \"q\", \"r\",\"t\", \"h\", \"c\"!")
    if(is.null(stepx))
      {
        if(method=="q"||method=="r"||method=="t"||method=="h")
          {
            stop("stepx is missing!")
          }
      }
    if(is.null(stepy))
      {
        if(method=="r")
          {
            stop("stepy is missing!")
          }
        else
          {
            stepy<-stepx
          }
      }
    if(is.null(gp))
      {
        if(method=="c")
          {
            stop("number of grid points is missing!")
          }
      }
    x.m<-max(x)-min(x)
    y.m<-max(y)-min(y)
    d.x<-max(x)-min(x)
    d.y<-max(y)-min(y)
    d.max<-max(d.x,d.y)
    if(method=="q"||method=="r")
      {     
        g.x<-min(x)-d.max*stepx
        g.y<-min(y)-d.max*stepy
        g1.x<-min(x)-d.max*stepx
        g1.y<-min(y)-d.max*stepy
        zx<-ceiling(d.x/d.max/stepx)+2
        for(i in 1:zx)
          {
            i<-i+1
            g1.x<-g1.x+d.max*stepx
            g.x<-c(g.x,g1.x) 
          }
         zy<-ceiling(d.y/d.max/stepy)+2
        for(i in 1:zy)
          {
            i<-i+1
            g1.y<-g1.y+d.max*stepy
            g.y<-c(g.y,g1.y) 
           }
        x.l<-length(g.x)
        y.l<-length(g.y)
        a1<-rep(g.x[1],y.l)
        for(i in 1:(x.l-1))
          {
            i<-i+1
            a1<-c(a1,rep(g.x[i],y.l))
          }
        a2<-rep(g.y,x.l)
        A<-matrix(c(a1,a2),nrow=x.l*y.l,ncol=2)
      }
    if(method=="t")
      {
        stepy<-sqrt(3)*stepx
        g.x<-min(x)-d.max*stepx
        g.y<-min(y)-d.max*stepy
        g1.x<-min(x)-d.max*stepx
        g1.y<-min(y)-d.max*stepy
        zx<-ceiling(d.x/d.max/stepx)+2
        for(i in 1:zx)
          {
            i<-i+1
            g1.x<-g1.x+d.max*stepx
            g.x<-c(g.x,g1.x) 
          }
        zy<-ceiling(d.y/d.max/stepy)+2
        for(i in 1:zy)
          {
            i<-i+1
            g1.y<-g1.y+d.max*stepy
            g.y<-c(g.y,g1.y) 
          }
        h.x<-min(x)-0.5*d.max*stepx
        h.y<-min(y)-0.5*d.max*stepy
        h1.x<-min(x)-0.5*d.max*stepx
        h1.y<-min(y)-0.5*d.max*stepy
        for(i in 1:(zx-1))
          {
            i<-i+1
            h1.x<-h1.x+d.max*stepx
            h.x<-c(h.x,h1.x) 
          }
        for(i in 1:(zy-1))
          {
            i<-i+1
            h1.y<-h1.y+d.max*stepy
            h.y<-c(h.y,h1.y) 
          }
        x.l1<-length(g.x)
        y.l1<-length(g.y)
        a11<-rep(g.x[1],y.l1)
        for(i in 1:(x.l1-1))
          {
            i<-i+1
            a11<-c(a11,rep(g.x[i],y.l1))
          }
        a21<-rep(g.y,x.l1)
        A1<-matrix(c(a11,a21),nrow=x.l1*y.l1,ncol=2)
        x.l2<-length(g.x)
        y.l2<-length(h.y)
        a12<-rep(h.x[1],y.l2)
        for(i in 1:(x.l2-1))
          {
            i<-i+1
            a12<-c(a12,rep(h.x[i],y.l2))
          }
         a22<-rep(h.y,x.l2)
        A2<-matrix(c(a12,a22),nrow=x.l2*y.l2,ncol=2)
        A<-rbind(A1,A2)
      }
    if(method=="h")
      {
        stepy<-2/3*sqrt(3)*stepx
        g.x<-min(x)-2/3*d.max*stepx
        g.y<-min(y)-d.max*stepy
        g1.x<-min(x)-2/3*d.max*stepx
        g1.y<-min(y)-d.max*stepy
        zx1<-ceiling((ceiling(d.x/d.max/stepx)+2)/2)
        for(i in 1:zx1)
          {
            i<-i+1
            g1.x<-g1.x+2*d.max*stepx
            g.x<-c(g.x,g1.x) 
          }
        zy1<-ceiling(d.y/d.max/stepy)+2
        for(i in 1:zy1)
          {
            i<-i+1
            g1.y<-g1.y+d.max*stepy
            g.y<-c(g.y,g1.y) 
          }
        x.l1<-length(g.x)
        y.l1<-length(g.y)
        a11<-rep(g.x[1],y.l1)
        for(i in 1:(x.l1-1))
          {
            i<-i+1
            a11<-c(a11,rep(g.x[i],y.l1))
          }
        a12<-rep(g.y,x.l1)
        A1<-matrix(c(a11,a12),nrow=x.l1*y.l1,ncol=2)         
        h.x<-min(x)
        h.y<-min(y)-d.max*stepy
        h1.x<-min(x)
        h1.y<-min(y)-d.max*stepy
        for(i in 1:(zx1-1))
          {
            i<-i+1
            h1.x<-h1.x+2*d.max*stepx
            h.x<-c(h.x,h1.x) 
          }
        for(i in 1:zy1)
          {
            i<-i+1
            h1.y<-h1.y+d.max*stepy
            h.y<-c(h.y,h1.y) 
          }
        x.l2<-length(h.x)
        y.l2<-length(h.y)
        a21<-rep(h.x[1],y.l2)
        for(i in 1:(x.l2-1))
          {
            i<-i+1
            a21<-c(a21,rep(h.x[i],y.l2))
          }
        a22<-rep(h.y,x.l2)
        A2<-matrix(c(a21,a22),nrow=x.l2*y.l2,ncol=2)
        k.x<-min(x)-d.max*stepx
        k.y<-min(y)-0.5*d.max*stepy
        k1.x<-min(x)-d.max*stepx
        k1.y<-min(y)-0.5*d.max*stepy
        for(i in 1:zx1)
          {
            i<-i+1
            k1.x<-k1.x+2*d.max*stepx
            k.x<-c(k.x,k1.x) 
          }
        for(i in 1:(zy1-1))
          {
            i<-i+1
            k1.y<-k1.y+d.max*stepy
            k.y<-c(k.y,k1.y) 
          }
        x.l3<-length(k.x)
        y.l3<-length(k.y)
        a31<-rep(k.x[1],y.l3)
        for(i in 1:(x.l3-1))
          {
            i<-i+1
            a31<-c(a31,rep(k.x[i],y.l3))
          }
        a32<-rep(k.y,x.l3)
        A3<-matrix(c(a31,a32),nrow=x.l3*y.l3,ncol=2)         
        l.x<-min(x)+1/3*d.max*stepx
        l.y<-min(y)-0.5*d.max*stepy
        l1.x<-min(x)+1/3*d.max*stepx
        l1.y<-min(y)-0.5*d.max*stepy
        for(i in 1:(zx1-1))
          {
            i<-i+1
            l1.x<-l1.x+2*d.max*stepx
            l.x<-c(l.x,l1.x) 
          }
        for(i in 1:(zy1-1))
          {
            i<-i+1
            l1.y<-l1.y+d.max*stepy
            l.y<-c(l.y,l1.y) 
          }
        x.l4<-length(l.x)
        y.l4<-length(l.y)
        a41<-rep(l.x[1],y.l4)
        for(i in 1:(x.l4-1))
          {
            i<-i+1
            a41<-c(a41,rep(l.x[i],y.l4))
          }
        a42<-rep(l.y,x.l4)
        A4<-matrix(c(a41,a42),nrow=x.l4*y.l4,ncol=2)
        A<-rbind(A1,A2,A3,A4)
      }
    if(method=="c")
      {
        g.x<-runif(gp,(min(x)-0.1*d.max),(max(x)+0.1*d.max))
        g.y<-runif(gp,(min(y)-0.1*d.max),(max(y)+0.1*d.max))
        A<-matrix(c(g.x,g.y),nrow=gp,ncol=2)
      }
    colnames(A)<-c("longitude","latitude")
    return(A)
  }
