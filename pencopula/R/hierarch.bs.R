 hierarch.bs <- function(x, d, plot.bsp, typ, penden.env, int=FALSE) #int optional! 
  {
    q <- get("q",penden.env)
    base <- get("base",penden.env)
    alpha <- get("alpha",penden.env)
    symmetric <- get("symmetric",penden.env)
    K <-  0
    ddb <- get("ddb",penden.env)#new
    
    knots <- knots.transform( d=K , alpha = get("alpha",penden.env),get("symmetric",penden.env))
  
    if(q==2) {help.dis <- knots[2]-knots[1]
      len.val <- length(knots)
      knots <- c(min(knots)-help.dis,knots,max(knots)+help.dis)
    }

    h <- 1
    if(base=="B-spline") {
      B.tilde <- my.bspline(h=h,q=q+1,knots=knots,y=x,K=length(knots),plot.bsp=plot.bsp,typ=typ)$base.den
      #if(q==2) B.tilde <- B.tilde[,-c(1,dim(B.tilde)[2])]
    }
    if(base=="Bernstein") {
      B.tilde <- apply(matrix(0:(ddb-1)),1,bernstein,x,n=(ddb-1))#[,get("k.order",penden.env)]
      dimBB <- dim(B.tilde)
    }
      
    #integriere B.tilde
    if(int& base=="B-spline") {
      if(q==1) index.h <- c(1,2)
      if(q==2) index.h <- c(1,2,3)
      int.B.tilde <- distr.func.help(B.tilde,knots,penden.env,q,y=x,index.h)
    }

    if(int&base=="Bernstein") {
      int.B.tilde <- int.bernstein(x,n=length(0:2^K))
    }

    if(base=="B-spline") {
      for ( K in 1:d)
        {
          h <- 1/(2**K)
          index <- (q-1) + 2*seq(1,2**(K-1),by=1)
          knots <- knots.transform(d=K, alpha = get("alpha",penden.env),get("symmetric",penden.env))

          if(q==2) {
            help.dis <- knots[2]-knots[1]
            len.val <- length(knots)
            knots <- c(min(knots)-help.dis,knots,max(knots)+help.dis)
          }
          
          BB <-  my.bspline(h,q=q+1,knots, y=x,K=length(knots),plot.bsp=plot.bsp,typ=typ)$base.den
              
          #integriere BB
        
          if(int & base=="B-spline") {
            index.h <- seq(max(index.h)+1,max(index.h)+length(index))
            BB.int <- distr.func.help(BB,knots,penden.env,q,y=x,index.h)
          }

        #if(int&base=="Bernstein") {
        #  BB.int <- int.bernstein(x,n=length(0:2^K))[,index]
        #}
        
          dimBB <- dim(BB)
       
          if(dimBB[1]>1) {
            B.tilde <-  cbind(B.tilde,BB[,index])
            if(int) int.B.tilde <- cbind(int.B.tilde,BB.int)
          }
          else {
            B.tilde <-  c(B.tilde,BB[,index])
            if(int) int.B.tilde <- c(int.B.tilde,BB.int)
          }
        }
      dimBB <- dim(BB)
      #jacobi <- matrix(1, length(x),1)#new
    }
   
    jacobi <- matrix(1, length(x),1)#new
    if(dimBB[1]>1) B.tilde.transform <- B.tilde * kronecker( matrix(jacobi) , matrix(1,1,dim(B.tilde)[2]))
    else B.tilde.transform <- B.tilde * kronecker( matrix(jacobi) , matrix(1,1,length(B.tilde)))
    if(int) {
      int.B.tilde.transform <- int.B.tilde * kronecker( matrix(jacobi) , matrix(1,1,dim(int.B.tilde)[2]))
      return(list(B.tilde=B.tilde.transform,int.B.tilde=int.B.tilde.transform))
    }
    else return(list(B.tilde=B.tilde.transform))
  }
    
