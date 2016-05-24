context("Graphics functions")

test_that("color", {
    cur <- palette()
    old <- lava:::mypal()
    expect_equivalent(col2rgb(cur),col2rgb(old))
    expect_equivalent(col2rgb(palette()),col2rgb(lava:::mypal(set=FALSE)))

    expect_equivalent(Col("red",0.5),rgb(1,0,0,0.5))
    expect_equivalent(Col(c("red","blue"),0.5),rgb(c(1,0),c(0,0),c(0,1),0.5))
    expect_equivalent(Col(c("red","blue"),c(0.2,0.5)),rgb(c(1,0),c(0,0),c(0,1),c(0.2,0.5)))
    expect_equivalent(Col(rgb(1,0,0),0.5),rgb(1,0,0,0.5))

    plot(0,xlim=c(0,1),ylim=c(0,1),type="n",ann=FALSE,axes=FALSE)
    devc1 <- devcoords()
    par(mar=c(0,0,0,0))
    plot(0,xlim=c(0,1),ylim=c(0,1),type="n",ann=FALSE,axes=FALSE)
    devc2 <- devcoords()
    figx <- c("fig.x1","fig.x2","fig.y1","fig.y2")
    devx <- c("dev.x1","dev.x2","dev.y1","dev.y2")
    expect_equivalent(devc1[figx],devc2[devx])    
    
})

if (requireNamespace("visualTest",quietly=TRUE) && requireNamespace("png",quietly=TRUE)) {

    gropen <- function(resolution=200,...) {
        tmpfile <- tempfile(fileext=".png")
        png(file=tmpfile,width=200,height=200)
        res <- dev.cur()
        return(structure(tmpfile,dev=res))
    }
    grcompare <- function(file1,file2,...) {
        res <- visualTest::isSimilar(file1,file2,...)
        unlink(c(file1,file2))
        return(res)
    }

    
    test_that("plotConf", {
        set.seed(1)
        x <- rnorm(50)
        y <- rnorm(50,x)
        z <- rbinom(50,1,0.5)
        d <- data.frame(y,z,x)
        l <- lm(y~x*z)
        
        d1 <- gropen()
        par(mar=c(0,0,0,0))
        plotConf(l,var2="z",col=c("black","blue"),alpha=0.5,legend=FALSE)
        dev.off()
        
        newd <- data.frame(x=seq(min(x),max(x),length.out=100))
        l0 <- lm(y~x,subset(d,z==0))
        ci0 <- predict(l0,newdata=newd,interval="confidence")
        l1 <- lm(y~x,subset(d,z==1))
        ci1 <- predict(l1,newdata=newd,interval="confidence")

        d2 <- gropen()
        par(mar=c(0,0,0,0))
        plot(y~x,col=c("black","blue")[z+1],pch=16,ylim=c(min(ci0,ci1,y),max(ci0,ci1,y)))
        lines(newd$x,ci0[,1],col="black",lwd=2)
        lines(newd$x,ci1[,1],col="blue",lwd=2)
        confband(newd$x,lower=ci0[,2],upper=ci0[,3],polygon=TRUE,col=Col("black",0.5),border=FALSE)
        confband(newd$x,lower=ci1[,2],upper=ci1[,3],polygon=TRUE,col=Col("blue",0.5),border=FALSE)
        points(y~x,col=c("black","blue")[z+1],pch=16)
        dev.off()
        
        expect_true(grcompare(d1,d2,threshold=5))

        
        d1 <- gropen()
        par(mar=c(0,0,0,0))
        l <- lm(y~z)
        plotConf(l,var2="z",var1=NULL,jitter=0,col="black",alpha=0.5,xlim=c(.5,2.5),ylim=range(y))
        dev.off()

        d2 <- gropen()
        par(mar=c(0,0,0,0))
        plot(y~I(z+1),ylim=range(y),xlim=c(0.5,2.5),pch=16,col=Col("black",0.5)) 
        l0 <- lm(y~-1+factor(z))
        confband(1:2,lower=confint(l0)[,1],upper=confint(l0)[,2],lwd=3,
                 center=coef(l0))
        dev.off()

        expect_true(grcompare(d1,d2,threshold=10))
                
  })


    test_that("forestplot", {
        set.seed(1)
        K <- 20
        est <- rnorm(K); est[c(3:4,10:12)] <- NA        
        se <- runif(K,0.2,0.4)        
        x <- cbind(est,est-2*se,est+2*se,runif(K,0.5,2))        
        rownames(x) <- unlist(lapply(letters[seq(K)],function(x) paste(rep(x,4),collapse="")))
        rownames(x)[which(is.na(est))] <- ""
        signif <- sign(x[,2])==sign(x[,3])        
        forestplot(x)
        ## TODO
    })

    test_that("plot.sim", {
        onerun2 <- function(a,b,...) {
            return(cbind(a=a,b=b,c=a-1,d=a+1))
        }
        R <- data.frame(a=1:2,b=3:4)
        val2 <- sim(onerun2,R=R,type=0)
        plot(val2)
        plot(val2,plot.type="single")
        density(val2)
        ## TODO
    })

    test_that("spaghetti", {
        K <- 5
        y <- "y"%++%seq(K)
        m <- lvm()
        regression(m,y=y,x=~u) <- 1
        regression(m,y=y,x=~s) <- seq(K)-1
        regression(m,y=y,x=~x) <- "b"
        d <- sim(m,5)
        dd <- mets::fast.reshape(d);
        dd$num <- dd$num+rnorm(nrow(dd),sd=0.5) ## Unbalance
        spaghetti(y~num,dd,id="id",lty=1,col=Col(1,.4),trend=TRUE,trend.col="darkblue")
        ## TODO
    })
    
    test_that("ksmooth", {
        ## TODO
    })

    test_that("plot.lvm", {
        ## TODO
        m <- lvm(y~u,u~x)
        latent(m) <- ~u
        plot(m)
        d <- sim(m,10)
        e <- estimate(m,d)
        plot(e)
        plot(lava:::beautify(m))
        g <- igraph.lvm(m)
        expect_true(inherits(g,"igraph"))
    })
    
    test_that("images", {
        ## TODO
    })

    test_that("labels,edgelabels", {
        ## TODO
    })

    test_that("colorbar", {
        ## TODO
    })

    test_that("fplot", {
        ## TODO
    })

    test_that("interactive", {
        ## TODO
    })

    test_that("pdfconvert", {
        ## TODO
    })

    test_that("plot.estimate", {
        ## TODO
    })

    test_that("logo", {
        lava(seed=1)
    })

    
    
    
}
