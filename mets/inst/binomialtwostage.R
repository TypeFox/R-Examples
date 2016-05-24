
onerunfam <- function(i,n,alr=0,manual=1,time=0,simplealr=1,theta=1) { ## {{{ 
### n=200; beta=0.2; theta=1; time=0; i=1
    print(i)
    dd <- simBinFam(n,beta=0,theta) 
    ddl <- fast.reshape(dd,varying="y",keep="y")
    out2t <- system.time(
        marg  <-  glm(y~+1,data=ddl,family=binomial())
        )
    ps <- predict(marg,type="response")
    if (time==1) print(out2t)

    if (manual==1) {
        if (time==1) print(date())
        ddl$ps <- ps
        fam <- familycluster.index(ddl$id)
        prtfam <- ddl[fam$familypairindex,]
        prtfam$subfam <- fam$subfamilyindex
### lave afhængighedsdesign pba af wide format zyg1*zyg2  feks
        prtfamclust <- data.frame(fast.reshape(prtfam,id="subfam"))
###     des <- model.matrix(~-1+factor(num1):factor(num2),data=prtfamclust)
        mf <- with(prtfamclust,(num1=="m")*(num2=="f")*1)
        mb <- with(prtfamclust,(num1=="m" | num1=="f")*(num2=="b1" | num2=="b2")*1)
        bb <- with(prtfamclust,(num1=="b1" )*(num2=="b1" | num2=="b2")*1)
        des <- cbind(mf,mb,bb)*1
        mulig <- (apply(des,2,sum)>0)
        names <- colnames(des)
        prtfamclust <- cbind(prtfamclust,des)
        prtfam <- fast.reshape(prtfamclust,varying=c("y","ps","num"),keep=c("y","ps","num","subfam","id",names))
        prtfam$famclust <- prtfam$id1
        destheta <- prtfam[,names]
        if (time==1) print(date())
        udt <-  system.time(
            udf <- binomial.twostage(prtfam$y,data=prtfam,
                                     clusters=prtfam$subfam, detail=0,
###	   score.method="nlminb",
                                     score.method="fisher.scoring",
                                     theta.des=prtfam[,names],
                                     max.clust=1000,iid=1,
                                     Nit=60,marginal.p=prtfam$ps,se.clusters=prtfam$famclust)
            )
        if (time==1) print(udt)

        zfam <- rbind(c(1,0,0), ## m-f
                      c(0,1,0),  ## m-b1
                      c(0,1,0),  ## m-b2
                      c(1,0,0), ## f-m
                      c(0,1,0),  ## f-b1
                      c(0,1,0),  ## f-b2
                      c(0,1,0), ## b1-m
                      c(0,1,0), ## b1-f
                      c(0,0,1), ## b1-b2
                      c(0,1,0), ## b2-m
                      c(0,1,0), ## b2-f
                      c(0,0,1)) ## b2-b1


        if (alr==1) {
	    if (simplealr==0) {
###       cvec <- (ddl$num=="m"| ddl$num=="f")*1 +(ddl$num=="b1"| ddl$num=="b2")*2
###       k <- 2
###       dmat <- rbind(c(1,0,0),c(0,1,0),c(0,0,1))
###       udz <- class2z(cvec,ddl$id,k,dmat)
###
###    ZMAST <- rep(1,12)
###    ZMAST <- cbind(ZMAST,c(0,0,0,0,1,1,0,1,1,0,1,1))
###
###   outl <- alr(ddl$y~+1,id=ddl$id,depmodel="general",ainit=rep(0.01,3),z=udz$z,zmast=0)
                if (!require(alr)) stop("'alr' package required")
                out4t <-  system.time(
                    outl <- alr(ddl$y~+1,id=ddl$id,depmodel="general",zlocs=rep(1:4,n),ainit=rep(0.01,3),z=zfam,zmast=1)
                    )
                if (time==1) print(out4t)

                outl <- c(summary(outl)$alpha[,1],summary(outl)$alpha[,2])
                names(outl) <- c(rep("alr",3),rep("se-alr",3))
	    } else {
                outl <- alr(ddl$y~+1,id=ddl$id,depmodel="exchangeable",ainit=rep(0.01,1))
                outl <- c(summary(outl)$alpha[,1],summary(outl)$alpha[,2])
                names(outl) <- c(rep("alr",1),rep("se-alr",1))
	    }

        }

    } else { ### med design formula
        form <- ~factor(num1)*factor(num2)
        udbin <- easy.binomial.twostage(marg,data=ddl,
                                        response="y",id="id",theta.formula=form,
                                        marginal.p=ps,
                                        score.method="fisher.scoring")

    }

###   if (alr==1) { ### alr til simpelt design
###   outl <- alr(ddl$y~ddl$x,id=ddl$id,depm="exchangeable", ainit=0.01)
###   outl <- summary(outl)$alpha
###   ud <- c(udbin$theta,udbin$var.theta^.5,udbin$hessi^.5,c(outl)[1:2])
###   names(ud) <- c("TWO","se-two","se-twoR","alr","se-alr")
###   } 
    ud <- c(udf$theta,diag(udf$var.theta)^.5)
    if (alr==1)  ud <- c(ud,outl)
    return(ud)
} ## }}} 

onerunfam2 <- function(i,n,alr=0,manual=1,time=0,theta=1) { ## {{{ 
### n=1000; beta=0.2; theta=1; time=0; i=1
    print(i)
    dd <- simBinFam(n,beta=0,theta) 
    ddl <- fast.reshape(dd,varying="y",keep="y")

    desfs <- function(x,num1="num1",num2="num2")
        { ## {{{ 
            mf <- (x[num1]=="m")*(x[num2]=="f")*1
            mb <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
            bb <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
            c(mf,mb,bb)
        } ## }}} 

    ud <- easy.binomial.twostage(y~+1,data=ddl,
                                 response="y",id="id",
                                 score.method="fisher.scoring",deshelp=0,
                                 theta.formula=desfs,desnames=c("pp","pc","cc"))

    ud <- c(ud$theta[,1],diag(ud$var.theta)^.5)

    zfam <- rbind(c(1,0,0), ## m-f
	      	  c(0,1,0),  ## m-b1
	      	  c(0,1,0),  ## m-b2
                  c(1,0,0), ## f-m
	      	  c(0,1,0),  ## f-b1
	      	  c(0,1,0),  ## f-b2
	      	  c(0,1,0), ## b1-m
	      	  c(0,1,0), ## b1-f
	      	  c(0,0,1), ## b1-b2
	      	  c(0,1,0), ## b2-m
	      	  c(0,1,0), ## b2-f
	      	  c(0,0,1)) ## b2-b1

    if (alr==1) {
        if (!require(alr)) stop("'alr' package required")
        outl <- alr(ddl$y~+1,id=ddl$id,
                    depmodel="general",zlocs=rep(1:4,n),ainit=rep(0.01,3),z=zfam,zmast=1)
        outl <- c(summary(outl)$alpha[,1],summary(outl)$alpha[,2])
        names(outl) <- c(rep("alr",3),rep("se-alr",3))
    }

    if (alr==1)  ud <- c(ud,outl)
    return(ud)
} ## }}} 
