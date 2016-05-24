tr <- function(y,tt,n){

   yf <- matrix(0,nrow=n,ncol=tt)
   for (i in 1:n) {
       yf[i,]<-y[(1+(i-1)*tt):(tt*i)]
   }
   yfm <- yf- colMeans(t(yf))
   yfm <- yfm[,1:(tt-1)]
   out <- matrix(t(yfm),nrow=nrow(yfm)*ncol(yfm),ncol=1)
   out
}


sse_calc <- function(y,x){
         #e <- y-x%*%qr.solve(x,y)
dat=data.frame(y,x)
myF=as.formula(paste(paste("y~", paste(names(dat[,-1]),collapse= "+")),"-1",sep=""))
         e<-resid(lm(myF,data=dat))
         out <- t(e)%*%e 
         out
} 

thr_sse <- function(y,q,r,cf,xt,ct,thresh,tt,n){
      nq <- nrow(q)
      sse <- matrix(0,nq,1) 

      for (qi in 1:nq){
            if (r[1]==0) {rr <- q[qi]
            } else { rr <- rbind(r,q[qi])}
            rr <- as.matrix(sort(rr))          
            xx <- cbind(xt,ct)
            tmp1=NULL; 

            for (j in 1:nrow(rr)){
                d <- (thresh < rr[j]);              
            for (i in 1:ncol(cf)) 
            {tmp1=cbind(tmp1,tr(cf[,i]*d,tt,n))}
             xx <- cbind(xx,tmp1)
            }
            sse[qi] <- sse_calc(y,xx)
        }
        sse
}

r_est <- function(y,r,trim,tt,qq1,qn1,qn,n,cf,xt,ct,thresh){
      if (max(r)==0){
          qq <- qq1;
          rr <- 0;
      } else {rr <- as.matrix(sort(r))
            i <- as.matrix(seq(1,qn1,by=1))
            nn <- colSums(qq1%*%matrix(1,1,nrow(rr))< matrix(1,nrow(qq1),1)%*%t(rr))
            nn <- as.matrix(nn)
            qnt <- qn*trim
            ii1 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn+qnt))
            ii2 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn-qnt))
            ii <- (ii1-ii2)%*%matrix(1,nrow(rr),1)
            qq <- as.matrix(qq1[ii!=1])
      }
      sse <- thr_sse(y,qq,rr,cf,xt,ct,thresh,tt,n)
      rihat <- which.min(sse)
      list(sse_b=sse[rihat],rhat_b=qq[rihat])
}

model <- function(r,trim,rep,it,qq1,cf,xt,ct,thresh,tt,qn1,n,qn,cc,yt,ty,k){
vgraph=1
      if (max(r)==0){
          qq <- qq1;
          rr <- 0;
      } else {rr <- as.matrix(sort(r))
            i <- as.matrix(seq(1,qn1,by=1))
            nn <- colSums(qq1%*%matrix(1,1,nrow(rr))< matrix(1,nrow(qq1),1)%*%t(rr))
            nn <- as.matrix(nn)
            qnt <- qn*trim
            ii1 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn+qnt))
            ii2 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn-qnt))
            ii <- (ii1-ii2)%*%matrix(1,nrow(rr),1)
            qq <- as.matrix(qq1[ii!=1])
      }
      sse <- thr_sse(yt,qq,rr,cf,xt,ct,thresh,tt,n)
      rihat <- which.min(sse)
      rhat <- qq[rihat]
      sse1 <- sse[rihat]
      lr <- (sse/sse1 - 1)*ty
      rhats <- as.matrix(qq[(lr < cc)])
      if (vgraph==1){
          if (it==0){
              titname=rbind("Figure 1","Confidence Interval Construction in Single Threshold Model")
              xname="Threshold Parameter"}
          if (it==1){
              titname=rbind("Figure 3","Confidence Interval Construction in Double Threshold Model")
              xname="First Threshold Parameter"}
          if (it==2){
              titname=rbind("Figure 2","Confidence Interval Construction in Double Threshold Model")
              xname="Second Threshold Parameter"}
          if (it==3){
              titname="Confidence Interval Construction in Triple Threshold Model"
              xname="Thrid Threshold Parameter"}
          yname="Likelihood Ratio"
          dev.new() 
          xxlim <- range(qq)
          yylim <- range(rbind(lr,cc))
          plot(qq,lr,lty=1,col=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
          lines(qq,matrix(1,nrow=nrow(qq),ncol=1)*cc,lty=2,col=2)
          title(main=titname,xlab=xname,ylab=yname)
      }
      if (max(r) != 0){
          cat ("Fixed Thresholds       ", t(rr), "\n")
          rrr <- sort(rbind(rr,rhat))
      } else { rrr <- rhat }
      rrr <- as.matrix(rrr)
      cat ("Threshold Estimate     ", rhat, "\n")
      cat ("Confidence Region      ", cbind(min(rhats),max(rhats)), "\n")
      cat ("Sum of Squared Errors  ", sse1, "\n")
      cat ("Trimming Percentage    ", trim, "\n")
      cat ("\n")
      cat ("\n")
      nr <- nrow(rrr)

      xx <- xt
      dd <- matrix(0,nrow=nrow(as.matrix(thresh)),ncol=nr)
      for (j in 1:nr){
          dd[,j] <- (thresh < rrr[j])
          d <- dd[,j]
          if (j>1) d = d - dd[,(j-1)];
          tmp2=NULL 
          for (i in 1:ncol(cf)) 
              {tmp2=cbind(tmp2,tr(cf[,i]*d,tt,n))};
          xx<- cbind(xx,tmp2) 
      }      
      
      d <- 1-dd[,nr]
      tmp3=NULL
      for (i in 1:ncol(cf))
      {tmp3=cbind(tmp3,tr(cf[,i]*d,tt,n))}
      
      xx<- cbind(xx,tmp3)

      xxi <- solve(t(xx)%*%xx)
      beta <- xxi%*%(t(xx)%*%yt)
      e <- yt - xx%*%beta
      xxe <- xx*(e%*%matrix((1),nrow=1,ncol=ncol(xx)))
      xxe <- t(xxe)%*%xxe
      sehet <- as.matrix(sqrt(diag(xxi%*%xxe%*%xxi)))
      sehomo <- as.matrix(sqrt(diag(xxi*as.vector((t(e)%*%e))/(ty-n-ncol(xx)))))
      beta <- cbind(beta,sehomo,sehet)
      cat ("Thresholds", "\n")
      cat (t(rrr), "\n")
      cat ("\n")
      cat ("Regime-independent Coefficients, standard errors, het standard errors,and t-stat", "\n")
      beta=cbind(beta,beta[,1]/beta[,3])
      cat(" Coeff",  "      std",  "        White","      tstat", "\n")
      beta <- format(beta, digits = 4, scientific = FALSE)
      for (j in 1:k) cat (beta[j,], "\n")
      cat ("\n")
      cat ("Regime-dependent Coefficients, standard errors, het standard errors,and t-stat", "\n")

      for (j in (k+1):nrow(beta)) cat (beta[j,], "\n")
      cat ("\n")
      cat ("\n")
      if (rep > 0){
          xx <- cbind(xt,ct)
          if (max(rr) != 0){
              for (j in 1:nrow(rr))
      tmp4=NULL;for (i in 1:ncol(cf)) {tmp4=cbind(tmp4,tr(cf[,i]*(thresh < rr[j]),tt,n))} ;xx <-cbind(xx,tmp4) 


          }
          yp <- xx%*%qr.solve(xx,yt)
          e <- yt-yp
          sse0 <- t(e)%*%e
          lrt <- (sse0/sse1-1)*ty
          cat ("LR Test for threshold effect  ", lrt, "\n")
          cat ("\n")
          cat ("\n")
          stats <- matrix(c(0),nrow=rep,ncol=1)
          for (j in 1:rep){
              eb <- matrix(c(0),nrow=n,ncol=(tt-1))
              for (i in 1:n) {
                  eb[i,]<-e[(1+(i-1)*(tt-1)):((tt-1)*i)]
              }
              yeb <- t(eb[ceiling(runif(n)*n),])
              yb <- yp + matrix(yeb,nrow=nrow(yeb)*ncol(yeb),ncol=1)
              sse0 <- sse_calc(yb,cbind(xt,ct))
              out <- r_est(yb,0,trim,tt,qq1,qn1,qn,n,cf,xt,ct,thresh)
              sse1 <- out$sse_b
              rhat_b <- out$rhat_b
              rrr <- rhat_b
              if (max(r) != 0){
                  for (jj in 1:length(r)){
                      sse0 <- sse1
                      out <- r_est(yb,rrr,trim,tt,qq1,qn1,qn,n,cf,xt,ct,thresh)
                      sse1 <- out$sse_b
                      rhat_b <- out$rhat_b
                      rrr <- rbind(rrr,rhat_b)
                  }
              }
              lrt_b <- (sse0/sse1-1)*ty
              stats[j] <- lrt_b
              cat ("Bootstrap Replication ", j, lrt_b, "\n")
          }
          cat ("\n")
          cat ("\n")
          stats <- as.matrix(sort(stats))
          crits <- as.matrix(stats[ceiling(rbind(.90,.95,.99)*rep)])
          cat ("Number of Bootstrap replications   ", rep, "\n")
          cat ("Bootstrap p-value                  ", mean(stats > as.vector(lrt)), "\n")
          cat ("Critical Values   ", crits[1], crits[2], crits[3], "\n")
          cat ("\n")
          cat ("\n")
      }
      rhat
}

ptm <- function(dep,ind1,ind2,d,bootn,trimn,qn,conf_lev,max_lag,t,n){

boot_1=bootn[1];boot_2=bootn[2];boot_3=bootn[3]
trim_1=trimn[1];trim_2=trimn[2];trim_3=trimn[3]

tt <- t-max_lag
ty <- n*(t-max_lag-1)

y  <- dep
cf=ind1
ct=NULL
for (j in 1:ncol(cf)) { ct=cbind(ct,tr(cf[,j],tt,n)) }


d1 <- d   # set to threshold variable
yt <- tr(y,tt,n)
#q1 <- ind2
#x <- cbind(q1,(q1^2),(q1^3),d1,(q1*d1))
x<-ind2
k <- ncol(x)
xt <- matrix(c(0),nrow=nrow(yt),ncol=k)
for (j in 1:k) { xt[,j]=tr(x[,j],tt,n) }
thresh <- d1
dd <- unique(thresh)
dd <- as.matrix(sort(dd))
qnt1 <- qn*trim_1
sq <- as.matrix(seq(trim_1,trim_1+(1/qn)*(qn-2*qnt1),by=1/qn))
qq1 <- as.matrix(dd[floor(sq*nrow(dd))])
qn1 <- nrow(qq1)
cc <- -2*log(1-sqrt(conf_lev))

#for (i in 1:1){
cat ("Number of Firms        ", n, "\n")
cat ("Number of Years used   ", tt, "\n")
cat ("Total Observations     ", ty, "\n")
cat ("Number of Quantiles    ", qn, "\n")
cat ("Confidence Level       ", conf_lev, "\n")
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")
cat ("Zero Threshold Model", "\n")
sse0 <- sse_calc(yt,cbind(xt,ct))
cat ("Sum of Squared Errors                   ", sse0, "\n")
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")

cat ("Single Threshold Model",  "\n")
cat ("\n")
rhat1 <- model(0,trim_1,boot_1,0,qq1,cf,xt,ct,thresh,tt,qn1,n,qn,cc,yt,ty,k)
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")

cat ("Double Threshold Model", "\n")
cat ("Trimming Percentage    ", trim_2, "\n")
cat ("\n")
cat ("First Iteration", "\n")
rhat2 <- model(rhat1,trim_2,boot_2,2,qq1,cf,xt,ct,thresh,tt,qn1,n,qn,cc,yt,ty,k)
cat ("Second Iteration", "\n")
rhat1 <- model(rhat2,trim_2,0,1,qq1,cf,xt,ct,thresh,tt,qn1,n,qn,cc,yt,ty,k)
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")

cat ("Triple Threshold Model", "\n")
cat ("Trimming Percentage    ", trim_3, "\n")
cat ("\n")
rhat3 <- model(rbind(rhat1,rhat2),trim_3,boot_3,3,qq1,cf,xt,ct,thresh,tt,qn1,n,qn,cc,yt,ty,k)
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")
#}
}


