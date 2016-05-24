                                        #
triangular.test.norm<-function(x, y=NULL, mu0=NULL, mu1, mu2=NULL,
                               delta=NULL, sigma=NULL, sigma2=NULL, 
                               alpha=0.05, beta=0.1,plot=TRUE){
  
  if(!is.null(mu2) && !is.null(mu0)){
    kind <- "two-sided"
  }else{
    kind <- "one-sided"
    ## FIXME:
    if(is.null(mu2))
      mu2 <- 0
    if(is.null(mu0))
      mu0 <- 0
  }
  ## initialize data structure
  
  if(is.null(y)){
    ## one sample
    n<-length(x)
    m<-NULL
    sample <- "one"
    if(!is.null(sigma)){
      variance="known"
    } else {
      variance="unknown"
      sigma <- sqrt(var(x))
    }
    # 
  } else {
    ## two sample
    sample <- "two"
    if(missing(mu1)){ #????
      mu1 <- mean(x)
    }
    if(missing(mu2)){
      if(!is.null(delta)){
        mu2 <- mu1+delta
      }else{
        mu2 <- mean(y)
      }
    }

    n<-length(x)
    m<-length(y)
    
    if(!is.null(sigma)){
      variance="known"
    } else {
      variance="unknown"
      if(n>=2)
        sigma1q<-var(x)
      else
        sigma1q<-Inf
      if(m>=2)
        sigma2q<-var(y)
      else
        sigma2q<-Inf
      sigma<-sqrt(((n-1)*sigma1q+(m-1)*sigma2q)/(n+m-2))
    }
    }
  if(sample=="one")
    a <- (1+qnorm(1-beta)/qnorm(1-alpha))*log(1/(2*alpha))/ ((mu1-mu0)/sigma)
  else
    a <- (1+qnorm(1-beta)/qnorm(1-alpha))*log(1/(2*alpha))/ ((mu1-mu2)/sigma)
  if(sample=="one")
    b <- ((mu1-mu0)/sigma)/(2*(1+qnorm(1-beta)/qnorm(1-alpha)))
  else
    b <- ((mu1-mu2)/sigma)/(2*(1+qnorm(1-beta)/qnorm(1-alpha)))
### cat("\n")
### cat(paste("a: ",a,"\n"))
### cat(paste("b: ",b,"\n"))
  A <- NULL
  B <- NULL
  if(kind=="two-sided"){
    ## additionally:
    if(sample=="one")
      A <- -(1+qnorm(1-beta)/qnorm(1-alpha))*log(1/(2*alpha))/ ((mu0-mu2)/sigma)
    else
      A <- (1+qnorm(1-beta)/qnorm(1-alpha))*log(1/(2*alpha))/ ((mu1-mu0)/sigma)
    if(sample=="one")
      B <- -((mu0-mu2)/sigma)/(2*(1+qnorm(1-beta)/qnorm(1-alpha)))
    else
      B <- ((mu1-mu0)/sigma)/(2*(1+qnorm(1-beta)/qnorm(1-alpha)))
### cat("\n")
### cat(paste("A: ",A,"\n"))
### cat(paste("B: ",B,"\n"))
  }

  ## build object:
  
  ret<-list(x=NULL,y=NULL,n=0,m=0,alpha=alpha,beta=beta,
            dist="normal", sample=sample, kind=kind,
            mu0=mu0, mu1=mu1, mu2=mu2,
            sigma=sigma, sigma2=sigma2, delta=delta,
            a=a,b=b,A=A,B=B,variance=variance,
            ## initially NULL:
            vn=NULL, zn=NULL, result=NULL, step=0)
  class(ret)<-"triangular.test"

  ## do initial update, this fills ret$x, ret$y and ret$zn and ret$vn:
  if(ret$sample=="one"){
    for(xi in x){
      ret <- update(ret,x=xi,initial=TRUE)
      if(ret$result!="continue")
        break
    }
  } else {
    l <- min(n,m)
    ret <- update(ret,x=x[1],y[1],initial=TRUE)
    if(l>1){
      for(i in 2:l){
        ret <- update(ret,x=x[i],initial=TRUE)
        if(ret$result!="continue")
          break
        ret <- update(ret,y=y[i],initial=TRUE)
        if(ret$result!="continue")
          break
      }
    }
    if(n<m && ret$result=="continue"){
      for(i in (l+1):m){
        ret <- update(ret,y=y[i],initial=TRUE)
        if(ret$result!="continue")
          break
      }
    }
    if(n>m && ret$result=="continue"){
      for(i in (l+1):n){
        ret <- update(ret,x=x[i],initial=TRUE)
        if(ret$result!="continue")
          break
      }
    }
  }

  if(plot){
    plot(ret)
    print(ret)
    
  }
  ret
}


update.triangular.test <- function(object, x=NULL, y=NULL, initial=FALSE, plot="last", recursive=FALSE, ...){
  if(!inherits(object,"triangular.test")){
    stop("object is no triangular test object")
  }
  if(length(x)>1 || length(y)>1){
    ## do recursive calls
    if(object$sample=="one"){
      for (xi in x){
        object <- update(object,xi,recursive=TRUE,plot=plot)
        if(object$result!="continue")
          break
      }
      ret <- object
      if(plot=="last"){
        plot(ret)
        print(ret)
      }
    } else {
      n <- length(x)
      m <- length(y)
      l <- min(n,m)

      if(l>0)
        for(i in 1:l){
          object <- update(object,x=x[i],recursive=TRUE,plot=plot)
          if(object$result!="continue")
            break
          object <- update(object,y=y[i],recursive=TRUE,plot=plot)
          if(object$result!="continue")
            break
        }
      
      if(n<m && object$result=="continue"){
        for(i in (l+1):m){
          object <- update(object,y=y[i],recursive=TRUE,plot=plot)
          if(object$result!="continue")
            break
        }
      }
      if(n>m && object$result=="continue"){
        for(i in (l+1):n){
          object <- update(object,x=x[i],recursive=TRUE,plot=plot)
          if(object$result!="continue")
            break
        }
      }
    }
    
    ret <- object
    if(plot=="last"){
      plot(ret)
      print(ret)
    }
  } else {
    ##single data given
    

    ## update data and sizes
    
    object$step <- object$step+1

    if(object$sample=="two"){
      if(is.null(x) && is.null(y))
        stop("both x an y are null!")
    } else {
      if(is.null(x))
        stop("x is null!")
    }
    
    if(!is.null(x)){
      object$x<-c(object$x,x)
      object$n<-object$n+1
    }
    if(object$sample=="two"){
      if(!is.null(y)){
        object$y<-c(object$y,y)
        object$m<-object$m+1
      } 
    }

    ## calculate a,b,vn,zn:

    if(object$dist=="normal"){
      ret <- update.triangular.test.norm(object,x,y,initial)
    }
                                        #  browser()
    if(object$dist=="bernoulli"){
      ret <- update.triangular.test.prop(object,x,y,initial)
    }

    ## use a,b,vn,zn for decision:
    ## default:
    ret$result="continue"
    
    if(is.nan(ret$a) | is.nan(ret$b)){
      ## continue if samples too small to estimate sigma and hence to
      ## calculate a and b:
      ret$result="continue"
    } else {
                                        #    if(!is.null(ret$zn) & !is.null(ret$vn)){
      if(is.na(ret$zn[ret$step]) | is.na(ret$vn[ret$step])){
        ret$result="continue"
                                        #      browser()
                                        #      }
      } else {
        ret$result <- decide.triangular.test(ret)
      }
    }
    
    
    if((!initial & !recursive) | plot=="all"){
      plot(ret)  
      print(ret)
    }
  } # end single / multiple data
  ret  
}

decide.triangular.test <- function(obj){
  if(obj$dist=="normal"){
    if((obj$mu0>obj$mu1 & obj$sample=="one")|
       (obj$mu2>obj$mu1 & obj$sample=="two")){
      parorder <- "ascending"
    }else{
      parorder <- "descending"
    }
  }
  if(obj$dist=="bernoulli"){
    if((obj$p0>obj$p1 & obj$sample=="one")|
       (obj$p2>obj$p1 & obj$sample=="two")){
      parorder <- "ascending"
    }else{
      parorder <- "descending"
    }

  }

  if(obj$kind=="one-sided"){
    ## one sided, means one triangle
      ## triangle points downward
    if(parorder=="ascending"){
      dir <- "down"
      upper <- -obj$a + 3*obj$b*obj$vn[obj$step]
      lower <-  obj$a +   obj$b*obj$vn[obj$step]
    } else {
      dir <- "up"
      upper <-  obj$a +   obj$b*obj$vn[obj$step]
      lower <- -obj$a + 3*obj$b*obj$vn[obj$step]
      
    }
    
    if(upper > obj$zn[obj$step] &
       obj$zn[obj$step] > lower)
      obj$result="continue"
    
    if((obj$zn[obj$step] <=  lower && dir=="down") |
       (obj$zn[obj$step] >=  upper && dir=="up"))
      obj$result="H1"
    
    if((obj$zn[obj$step] <=  lower && dir=="up") |
       (obj$zn[obj$step] >=  upper && dir=="down"))
      obj$result="H0"
  } else {
    ## two sided
    ## means two triangles

    if(parorder=="ascending"){
      ## triangle1 points downward
      dir <- "down"
      upper <- -obj$a + 3*obj$b*obj$vn[obj$step]
      lower <-  obj$a +   obj$b*obj$vn[obj$step]
      Dir <- "up"
      Upper <-  obj$A +   obj$B*obj$vn[obj$step]
      Lower <- -obj$A + 3*obj$B*obj$vn[obj$step]
    } else {
      dir <- "up"
      upper <-  obj$a +   obj$b*obj$vn[obj$step]
      lower <- -obj$a + 3*obj$b*obj$vn[obj$step]
      Dir <- "down"
      Upper <- -obj$A + 3*obj$B*obj$vn[obj$step]
      Lower <-  obj$A +   obj$B*obj$vn[obj$step]
      
    }
    
                                        
    if((upper > obj$zn[obj$step] &
        obj$zn[obj$step] > lower) |
       (Upper > obj$zn[obj$step] &
        obj$zn[obj$step] > Lower)){
      obj$result="continue"
    }else{   
      if(((obj$zn[obj$step] <=  lower && dir=="down") |
          (obj$zn[obj$step] >=  upper && dir=="up")) |
         ((obj$zn[obj$step] <=  Lower && Dir=="down") |
          (obj$zn[obj$step] >=  Upper && Dir=="up")))
        obj$result="H1"
      
      if(((obj$zn[obj$step] <=  lower && dir=="up") |
          (obj$zn[obj$step] >=  upper && dir=="down")) &
         ((obj$zn[obj$step] <=  Lower && Dir=="up") |
          (obj$zn[obj$step] >=  Upper && Dir=="down")))
        obj$result="H0"
    }
  }
obj$result
}



update.triangular.test.norm <- function(object, x=NULL, y=NULL, initial=FALSE, plot="last", recursive=FALSE, ...){  

  if(!inherits(object,"triangular.test"))
    stop("object is not a triangular test object!")
  if(!initial)
    if(object$result!="continue")
      stop("Triangular test was already finished!")
  
  ## updates:
  if(object$variance=="unknown" && !initial){
    ## update unknown variance:
    if(object$sample=="one"){
      ## one sample
      object$sigma <- sqrt(var(object$x))
      ## update a and b:
      object$a <- (1+qnorm(1-object$beta)/qnorm(1-object$alpha))*
        log(1/(2*object$alpha))/((object$mu0-object$mu1)/object$sigma)
      object$b <- ((object$mu0-object$mu1)/object$sigma)/
        (2*(1+qnorm(1-object$beta)/qnorm(1-object$alpha)))
      if(object$kind=="two-sided"){
        object$A <- (1+qnorm(1-object$beta)/qnorm(1-object$alpha))*
          log(1/(2*object$alpha))/((object$mu0-object$mu1)/object$sigma)
        object$B <- ((object$mu0-object$mu1)/object$sigma)/
          (2*(1+qnorm(1-object$beta)/qnorm(1-object$alpha)))
      }
    } else {
      ## two sample
      if(object$n>=2)
        sigma1q<-var(object$x)
      else
        sigma1q<-Inf
      if(object$m>=2)
        sigma2q<-var(object$y)
      else
        sigma2q<-Inf
      
      ## object$sigma<-sqrt(((object$n-1)*sigma1q+(object$m-1)*sigma2q)/(object$n+object$m-2))
      object$sigma <- sqrt((sum(object$x^2)+sum(object$y^2)-
                         ((sum(object$x)+sum(object$y))^2/(object$n+object$m)))/
                        (object$n+object$m))
      ## update a and b:
      object$a <- (1+qnorm(1-object$beta)/qnorm(1-object$alpha))*
        log(1/(2*object$alpha))/((object$mu1-object$mu2)/object$sigma)
      object$b <- ((object$mu1-object$mu2)/object$sigma)/
        (2*(1+qnorm(1-object$beta)/qnorm(1-object$alpha)))
      if(object$kind=="two-sided"){
        object$A <- (1+qnorm(1-object$beta)/qnorm(1-object$alpha))*
          log(1/(2*object$alpha))/((object$mu0-object$mu1)/object$sigma)
        object$B <- ((object$mu0-object$mu1)/object$sigma)/
          (2*(1+qnorm(1-object$beta)/qnorm(1-object$alpha)))
      }
    }
  }

  ## calculate zn, vn:
  
  if(object$sample=="one"){
    ## one sample
    object$zn[object$step] <- sum(object$x-object$mu0)/sqrt(sum((object$x-object$mu0)^2)/object$n)
    object$vn[object$step] <- object$n-object$zn[object$step]^2/(2*object$n)
  } else {
    ## two sample, example 5.10, sigma unknown?
    ## mean ? mu0/mu1 ?????
    Snsqrd <- (sum(object$x^2) + sum(object$y^2) -
               ((sum(object$x) + sum(object$y))^2/(object$n + object$m)))/(object$n + object$m)

    if(!is.na(Snsqrd) & Snsqrd!=0){
      object$zn[object$step] <- object$n*object$m/(object$n+object$m)*
        (mean(object$x)-mean(object$y))/sqrt(Snsqrd)
      object$vn[object$step] <- object$n*object$m/(object$n+object$m)-object$zn[object$step]^2/
        (2*(object$n+object$m))
    } else {
      object$zn[object$step] <- NA
      object$vn[object$step] <- NA
    }
  }
  
  
  object
}


plot.triangular.test <- function(x, ...){
  if(!inherits(x,"triangular.test"))
    stop("x is not a triangular test object!")
  
  if(is.nan(x$a) | is.nan(x$b) | is.infinite(x$a) | is.infinite(x$b) | all(is.na(x$vn))){
    warning("cannot plot, a, b unknown and not yet estimated!")
  } else {
    if(x$dist=="normal"){
      if((x$mu0>x$mu1 & x$sample=="one")|
         (x$mu2>x$mu1 & x$sample=="two")){
        parorder <- "ascending"
      }else{
        parorder <- "descending"
      }
    }
    if(x$dist=="bernoulli"){
      if((x$p0>x$p1 & x$sample=="one")|
         (x$p2>x$p1 & x$sample=="two")){
        parorder <- "ascending"
      }else{
        parorder <- "descending"
      }
    }

    ## if plot is possible:
    if(x$kind=="one-sided"){
      ## one sided, one triangle
      vright <- x$a/x$b
      ##      vmin <- if(!is.na(x$vn[1])) x$vn[1] else x$vn[2]
      vmin <- 0
      vmax <- vright
      if(x$sample=="one"){ # 
        ## one sample
        if(parorder=="ascending"){
          z1 <- -x$a+  x$b*vmin
          z2 <-  x$a+3*x$b*vmin
          zright <- 2*x$a
        } else {
          z1 <- -x$a+3*x$b*vmin 
          z2 <-  x$a+  x$b*vmin
          zright <-  2*x$a
        }
      } else {
        ## two sample
        if(parorder=="ascending"){
          z1 <- -x$a+3*x$b*vmin 
          z2 <-  x$a+  x$b*vmin
          zright <-  2*x$a
        } else {
          z1 <- -x$a+  x$b*vmin 
          z2 <-  x$a+3*x$b*vmin
          zright <- 2*x$a
        }
      }
    } else {
      ## two sided
### TODO
      ## two triangles
        ## one sample
      vright <- x$a/x$b
      Vright <- x$A/x$B
      ##vmin <- if(!is.na(x$vn[1])) x$vn[1] else x$vn[2]
      vmin <- 0
      vmax <- max(vright,Vright)
      
      if(x$sample=="one"){ # 
        if(parorder=="ascending"){
          z1 <- -x$a+  x$b*vmin
          z2 <-  x$a+3*x$b*vmin
          zright <- 2*x$a
          Z1 <- -x$A+3*x$B*vmin 
          Z2 <-  x$A+  x$B*vmin
          Zright <-  2*x$A
        } else {
          z1 <- -x$a+3*x$b*vmin 
          z2 <-  x$a+  x$b*vmin
          zright <-  2*x$a
          Z1 <- -x$A+  x$B*vmin
          Z2 <-  x$A+3*x$B*vmin
          Zright <- 2*x$A
        }
      } else {
        ## two sample
        if(parorder=="ascending"){
          z1 <- -x$a+3*x$b*vmin 
          z2 <-  x$a+  x$b*vmin
          zright <-  2*x$a
          Z1 <- -x$A+3*x$B*vmin 
          Z2 <-  x$A+  x$B*vmin
          Zright <-  2*x$A
        } else {
          z1 <- -x$a+  x$b*vmin 
          z2 <-  x$a+3*x$b*vmin
          zright <- 2*x$a
          Z1 <- -x$A+  x$B*vmin
          Z2 <-  x$A+3*x$B*vmin
          Zright <- 2*x$A
        }
      }
    }
    
    
    ## determine size of plot
    if(x$kind=="one-sided"){
      zmin <- min(c(z1,z2,zright),na.rm=TRUE)
      zmax <- max(c(z1,z2,zright),na.rm=TRUE)
    } else {
      zmin <- min(c(z1,z2,zright,Zright),na.rm=TRUE)
      zmax <- max(c(z1,z2,zright,Zright),na.rm=TRUE)
    }
    
    plot(0,0,xlim=c(vmin,vmax),ylim=c(zmin,zmax),type="n",
         main="Triangular Test", xlab="v_n",ylab="z_n")
    
    
    if(x$kind=="one-sided"){
      ## one sided, one triangle:
      polygon(c(vmin,vmin,vright,vmin),
              c(z1,z2,zright,z1),col="grey",border="black")
      ## labels, depending on relation between m0 and mu1 (mu1 /mu2)

      if(parorder=="ascending"){
        text(vmax, zmax, "H0")
        text(vmin, zmin, "H1")
      } else {
        text(vmin, zmax, "H1")
        text(vmax, zmin, "H0")
      }
    } else {
      ## two sided, two triangles:
      ## need to calculate interection of triangle sides for unsymmetric case:
      ##
      ##  z1
      ##    ' .
      ##  Z2.... d3 ..........Zright
      ##           .     . '   
      ##  z2.       .d1'
      ##      'd2.'   * .
      ##  Z1.'       '  .'
      ##                   ' zright
      ##  vmin         vright  Vright
      ##
      ## d1: intersection Z1-Zright x z2-zright
      ## d2: intersection Z1-Zright x z1-zright
      ## d3: intersection Z2-Zright x z2-zright

      d1x <- ((vright - vmin) * Vright * Z1 + (vmin^2  - vmin * vright) * Zright + (vmin *  Vright - vmin^2 ) * zright + (vmin * vright - vright * Vright) * z1)/((vright - vmin) * Z1 + (vmin - vright) * Zright + (Vright - vmin) * zright + (vmin - Vright) * z1)
      d1y <- (((Vright - vmin) * zright + (vright - Vright) * z1) * Z1 + (vmin - vright) * z1 * Zright)/((vright - vmin) * Z1 + (vmin - vright) * Zright + (Vright - vmin) * zright + (vmin - Vright) * z1)

      d2x <- ((vright - vmin) * Vright * Z1 + (vmin^2  - vmin * vright) * Zright + (vmin * Vright - vmin^2 ) * zright + (vmin * vright - vright * Vright) * z2)/((vright - vmin) * Z1 + (vmin - vright) * Zright + (Vright - vmin) * zright + (vmin - Vright) * z2)
      d2y <- (((Vright - vmin) * zright + (vright - Vright) * z2) * Z1 + (vmin - vright) * z2 * Zright)/((vright - vmin) * Z1 + (vmin - vright) * Zright + (Vright - vmin) * zright + (vmin - Vright) * z2)

      ## dont plot left of vmin
      if(d2x<vmin){
        d2x <- vmin
        d2y <- z2
      }
      
      d3x <- ((vright - vmin) * Vright * Z2 + (vmin^2  - vmin * vright) * Zright + (vmin * Vright - vmin^2 ) * zright + (vmin * vright - vright * Vright) * z1)/((vright - vmin) * Z2 + (vmin - vright) * Zright + (Vright - vmin) * zright + (vmin - Vright) * z1)
      d3y <- (((Vright - vmin) * zright + (vright - Vright) * z1) * Z2 + (vmin - vright) * z1 * Zright)/((vright - vmin) * Z2 + (vmin - vright) * Zright + (Vright - vmin) * zright + (vmin - Vright) * z1)

      ## dont plot left of vmin
      if(d3x<vmin){
        d3x <- vmin
        d3y <- Z2
      }
      polygon(c(vmin,vmin,d2x,vright,d1x,Vright,d3x,vmin),
              c(z1,Z1,d2y,zright,d1y,Zright,d3y,z1),
              col="grey",
              border="black")
              

      
      ## labels, depending on relation between m0 and mu1 (mu1 /mu2)

    
      text(vmin, zmax, "H1")
      text(vmin, zmin, "H1")
      text(vmax, (zmin+zmax)/2, "H0")
    }

    
    ## the line:
    lines(x$vn, x$zn)
    if(length(x$vn)<=20)
      points(x$vn, x$zn)
  } 
}

print.triangular.test <- function(x, ...){
  if(!inherits(x,"triangular.test"))
    stop("x is not a triangular test object!")
  cat(paste("Triangular Test for",x$dist,"distribution\n\n"))
  if(x$dist=="normal"){
    if(x$variance=="known")
      cat(paste("Sigma known:", x$sigma,"\n"))
    else
      cat(paste("Sigma unknown, estimated as ", x$sigma, "\n"))
    if(x$sample=="one"){
      ## one sample
      if(x$kind=="one-sided"){
        if(x$mu1>x$mu0){
          cat(paste("\nH0: mu=",x$mu0," versus H1: mu>",x$mu1,"\n"))
        }else{
          cat(paste("\nH0: mu=",x$mu0," versus H1: mu<",x$mu1,"\n"))
        }
      }else{
        if(x$mu2>x$mu1){
          cat(paste("\nH0: mu=",x$mu0," versus H1: mu>",x$mu2," or mu<",x$mu1,"\n"))
        }else{
          cat(paste("\nH0: mu=",x$mu0," versus H1: mu>",x$mu1," or mu<",x$mu2,"\n"))
        }

      }
        
    } else { # two sample
      if(x$kind=="one-sided"){
        if(x$mu2>x$mu1){
          cat(paste("\nH0: mu1=mu2=",x$mu1," versus H1: mu1=",x$mu1," mu2>=",x$mu2,"\n"))
        }else{
          cat(paste("\nH0: mu1=mu2=",x$mu1," versus H1: mu1=",x$mu1," mu2<=",x$mu2,"\n"))
        }
      } else {
        ## two sided
        if(x$mu2>x$mu1){
          cat(paste("\nH0: mu1=mu2=",x$mu1," versus H1: mu1=",x$mu1," and mu2>=",x$mu2," or mu2<=",x$mu0,"\n"))          
        } else {
          cat(paste("\nH0: mu1=mu2=",x$mu1," versus H1: mu1=",x$mu1," and mu2>=",x$mu0," or mu2<=",x$mu2,"\n"))          
        }
      }
    }
  }
  
  if(x$dist=="bernoulli"){

    if(x$sample=="one"){
      ## one sample
      if(x$kind=="one-sided"){
        if(x$p1>x$p0){
          cat(paste("H0: p=",x$p0," versus H1: p>",x$p1,"\n"))
        }else{
          cat(paste("H0: p=",x$p0," versus H1: p<",x$p1,"\n"))
        }
      }else{
        if(x$p2>x$p1){
          cat(paste("H0: p=",x$p0," versus H1: p>",x$p2," or p<",x$p1,"\n"))
        }else{
          cat(paste("H0: p=",x$p0," versus H1: p>",x$p1," or p<",x$p2,"\n"))
        }

      }
        
    }else{ # two sample
      if(x$kind=="one-sided"){
        if(x$p2>x$p1){
          cat(paste("H0: p1=p2=",x$p1," versus H1: p1=",x$p1," p2>=",x$p2,"\n"))
        }else{
          cat(paste("H0: p1=p2=",x$p1," versus H1: p1=",x$p1," p2<=",x$p2,"\n"))
        }
      }else{
        ## two sided
        if(x$p2>x$p1){
          cat(paste("H0: p1=p2=",x$p1," versus H1: p1=",x$p1," and p2>=",x$p2," or p2<=",x$p0,"\n"))          
        }else{
          cat(paste("H0: p1=p2=",x$p1," versus H1: p1=",x$p1," and p2>=",x$p0," or p2<=",x$p2,"\n"))          
        }
      }
    }
  }


  cat(paste("alpha:", x$alpha," beta:", x$beta,"\n\n"))
  if(x$result=="continue"){
    cat("Test not finished, continue by adding single data via update()\n")
    if(x$sample=="one"){
      cat(paste("current sample size for x: ",x$n,"\n"))
    } else {
      cat(paste("current sample size for x: ",x$n,"\n"))
      cat(paste("current sample size for y: ",x$m,"\n"))
    }

###   cat("zn:\n")
###   print(x$zn)
###   cat("vn:\n")
###   print(x$vn)
###    cat("\nBounds: al=",x$al,"< zn-b*vn =", x$zn-x$b*x$vn, "<", x$au,"=au\n")
  } else {
    cat(paste("Test finished: accept", x$result,"\n"))
    if(x$sample=="one"){
      cat(paste("Sample size for x: ",x$n,"\n"))
    } else {
      cat(paste("Sample size for x: ",x$n,"\n"))
      cat(paste("Sample size for y: ",x$m,"\n"))
    }
###   cat("zn:\n")
###   print(x$zn)
###   cat("vn:\n")
###   print(x$vn)
###       cat("\nBounds: zn-b*vn =", x$zn-x$b*x$vn," not in  [",x$al,",", x$au,"= [al,au]\n")
  }
}



triangular.test.prop<-function(x, y=NULL, p0=NULL, p1=NULL, p2=NULL, alpha=0.05, beta=0.1, delta=NULL,plot=TRUE){

  
  if(!is.null(p0) && !is.null(p2)){
    kind <- "two-sided"
  }else{
    kind <- "one-sided"
    # FIXME
    if(is.null(p2))
       p2 <- 0
    if(is.null(p0))
       p0 <- 0
  }


     
  
  if(is.null(y)){
    ## one sample
    n<-length(x)
    m<-NULL
    
    sample <- "one"
    if(any(x!=0 & x!=1))
      stop("only 0 or 1 allowed for x!")
    if(is.null(p1)){
      if(!is.null(delta))
        p1 <- p0+delta
      else
        stop("p1 missing and no delta given!")
      }
  } else {
    ## two sample
    sample <- "two"
    n<-length(x)
    m<-length(y)
    if(any(x!=0 & x!=1))
      stop("only 0 or 1 allowed for x!")
    if(any(y!=0 & y!=1))
      stop("only 0 or 1 allowed for y!")
    if(is.null(p0) && is.null(p2) && is.null(p1)){
      ## estimate ps
      prop <- "estimate"
      p1 <- sum(x==1)/length(x)
      p2 <- sum(y==1)/length(y)
    } else {
      prop <- "given"
    }
  }
  
  if(sample=="one"){
    theta1 <- log((p1*(1-p0))/(p0*(1-p1)))
  }else{
    theta1 <- log((p1*(1-p2))/(p2*(1-p1)))
  }
  a <- (1+qnorm(1-beta)/qnorm(1-alpha))*log(1/(2*alpha))/theta1
  b <- theta1/(2*(1+qnorm(1-beta)/qnorm(1-alpha)))
  A <- NULL
  B <- NULL
  if(kind=="two-sided"){
    ## additionally:
    if(sample=="one"){
      Theta1 <- log((p2*(1-p0))/(p0*(1-p2)))
    }else{
      Theta1 <- log((p1*(1-p0))/(p0*(1-p1)))
    }
    A <- (1+qnorm(1-beta)/qnorm(1-alpha))*log(1/(2*alpha))/Theta1
    B <- Theta1/(2*(1+qnorm(1-beta)/qnorm(1-alpha)))
  }
### TODO 
  ## build object:
  
  ret<-list(x=NULL,y=NULL,n=0,m=0,alpha=alpha,beta=beta,
            dist="bernoulli", sample=sample, kind=kind,
            p0=p0,p1=p1,p2=p2,
            a=a,b=b,A=A,B=B,
            ## initially NULL:
            vn=NULL, zn=NULL, result=NULL, step=0)
  class(ret)<-"triangular.test"
  
  ## do initial update, this fills ret$x, ret$y and ret$zn and ret$vn:
  if(ret$sample=="one"){
    for(xi in x){
      ret <- update(ret,x=xi,initial=TRUE)
      if(ret$result!="continue")
        break
    }
  } else {
    ret <- update(ret,x=x[1],y=y[1],initial=TRUE)
#      if(ret$result!="continue")
#        break
# FIXME?
    l <- min(n,m)
    if(l>1){
      for(i in 2:l){
        ret <- update(ret,x=x[i],initial=TRUE)
        if(ret$result!="continue")
          break
        ret <- update(ret,y=y[i],initial=TRUE)
        if(ret$result!="continue")
          break
      }
    }
    if(n<m & ret$result=="continue"){
      for(i in (l+1):m){
        ret <- update(ret,y=y[i],initial=TRUE)
        if(ret$result!="continue")
          break
      }
    }
    if(n>m & ret$result=="continue"){
      for(i in (l+1):n){
        ret <- update(ret,x=x[i],initial=TRUE)
        if(ret$result!="continue")
          break
      }
    }
  }
  if(plot) {
    plot(ret)
    print(ret)
  }
  ret
}

update.triangular.test.prop <- function(object, x=NULL, y=NULL, initial=FALSE, plot="last", recursive=FALSE, ...){
  if(!inherits(object,"triangular.test"))
    stop("object is not a triangular test object!")
  if(!initial)
    if(object$result!="continue")
      stop("Triangular test was already finished!")
  if(!is.null(x))
    if(any(x!=0)&any(x!=1))
      stop("only 0 or 1 for x!")
  if(!is.null(y))
    if(any(y!=0)&any(y!=1))
      stop("only 0 or 1 for y!")
 
  ## calculate zn, vn:
  
  if(object$sample=="one"){
    object$r <- sum(object$x)
    object$zn[object$step] <- object$r-object$n*object$p0
    object$vn[object$step] <- object$n*object$p0*(1-object$p0)
  } else {
    object$r1 <- sum(object$x)
    object$r2 <- sum(object$y)
    object$zn[object$step] <- (object$m*object$r2-object$n*object$r1)/(object$n+object$m)
    object$vn[object$step] <- object$n*object$m*(object$r1+object$r2)*
      (object$n+object$m-(object$r1+object$r2))/(object$n+object$m)^3
  }
  object
}
