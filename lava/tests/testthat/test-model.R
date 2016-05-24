context("Model specification")

test_that("Basic model building blocks", {
    m <- lvm(y[m]~x)
    covariance(m) <- y~z
    expect_true(covariance(m)$rel["z","y"]==1)
    expect_true(regression(m)$rel["x","y"]==1)

    ## Children parent,nodes
    expect_match(children(m,~x),"y")
    expect_match(parents(m,~y),"x")    
    expect_equivalent(parents(m),vars(m))
    expect_equivalent(children(m),vars(m))

    ## Remove association
    cancel(m) <- y~z+x
    expect_true(covariance(m)$rel["z","y"]==0)
    expect_true(regression(m)$rel["x","y"]==0)

    ## Remove variable
    kill(m) <- ~x
    expect_equivalent(vars(m),c("y","z"))
    expect_true(intercept(m)["y"]=="m")  

    m <- lvm(c(y1,y2,y3)~x)
    d <- sim(m,50)
    e <- estimate(m,d)
    ## Equivalence
    ##equivalence(e,silent=TRUE)


    ## formula
    f <- formula(m,all=TRUE)
    expect_true(length(f)==length(vars(m)))
    expect_true(all(unlist(lapply(f,function(x) inherits(x,"formula")))))

    ## Parametrization
    m <- lvm(c(y1,y2,y3)~u)
    latent(m) <- ~u
    m2 <- fixsome(m,param=NULL)
    expect_true(all(is.na(regression(m2)$values)))
    m2 <- fixsome(m,param="relative")
    expect_true(regression(m2)$values["u","y1"]==1)
    expect_true(intercept(m2)[["y1"]]==0)
    m2 <- fixsome(m,param="hybrid")
    expect_true(regression(m2)$values["u","y1"]==1)
    expect_true(intercept(m2)[["u"]]==0)
    m2 <- fixsome(m,param="absolute")
    expect_true(all(is.na(regression(m2)$values)))
    expect_true(intercept(m2)[["u"]]==0)
    expect_true(covariance(m2)$values["u","u"]==1)

    ## Merge
    m1 <- lvm(c(y1,y2,y3)~1*u1[m1:v1])
    latent(m1) <- ~u1
    m2 <- lvm(c(y1,y2,y3)~2*u2[m2:v2])
    latent(m2) <- ~u2
    mm <- m1%++%m2

    expect_true(covariance(mm)$labels["u1","u1"]=="v1")
    expect_true(intercept(mm)[["u2"]]=="m2")

    ## LISREL
    mm <- fixsome(mm)
    L <- lisrel(mm,rep(1,length(coef(mm))))
    expect_equivalent(L$B,matrix(0,2,2))
    expect_equivalent(L$Theta,diag(3))
    expect_equivalent(L$Psi,diag(2))
    
}) 


test_that("Linear constraints", {
    m <- lvm(c(y[m:v]~b*x))
    constrain(m,b~a) <- base::identity
})


test_that("Graph attributes", {
    require("graph")
    m <- lvm(y~x)
    g1 <- graph::updateGraph(plot(m,noplot=TRUE))
    m1 <- graph2lvm(g1)
    expect_equivalent(m1$M,m$M)
    
    col <- "blue"; v <- "y"
    g1 <- lava::addattr(g1,"fill",v,col)
    expect_match(col,graph::nodeRenderInfo(g1)$fill[v])
    nodecolor(m,v) <- "blue"
    g2 <- Graph(m,add=TRUE)
    expect_true(inherits(g2,"graph"))
    expect_match(col,graph::nodeRenderInfo(g2)$fill[v])
    expect_match(addattr(g2,"fill")["y"],"blue")
    graph::graphRenderInfo(g2)$rankdir <- "LR"
    Graph(m) <- g2
    expect_true(graph::graphRenderInfo(Graph(m))$rankdir=="LR")

    ## Labels
    labels(m) <- c(y="Y")
    addattr(Graph(m,add=TRUE),"label")
    expect_true(addattr(finalize(m),"label")["y"]=="Y")
    labels(g2) <- c(y="Y")
    expect_true(graph::nodeRenderInfo(g2)$label["y"]=="Y")

    edgelabels(m,y~x) <- "a"
    expect_true(!is.null(edgelabels(finalize(m))))
})

 
test_that("Categorical variables", {
    m <- lvm()
    categorical(m,K=3,p=c(0.1,0.5)) <- ~x
    d1 <- simulate(m,10,seed=1)
    categorical(m,K=3) <- ~x
    d2 <- simulate(m,10,seed=1)
    expect_false(identical(d1,d2))
    
    regression(m,additive=FALSE,y~x) <- c(0,-5,5)
    d <- simulate(m,100,seed=1)
    l <- lm(y~factor(x),d)
    expect_true(sign(coef(l))[2]==-sign(coef(l))[3])
     
})
