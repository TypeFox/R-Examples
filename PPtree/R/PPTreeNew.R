 PP.Tree<- function (PPmethod, i.class, i.data, weight = TRUE, r = NULL,  lambda = NULL, cooling = 0.999, temp = 1, energy = 0.01,  ...) 
{
    i.data <- as.matrix(i.data)
    Find.proj <- function(i.class, i.data, PPmethod, r, lambda,  ...) 
   {
        n <- nrow(i.data)
        p <- ncol(i.data)
        g <- table(i.class)
        g.name <- as.numeric(names(g))
        G <- length(g)
        a <- PP.optimize.anneal(PPmethod, 1, i.data, i.class, 
            std = TRUE, cooling, temp, energy, r, lambda)
        proj.data <- as.matrix(i.data) %*% a$proj.best
        sign <- sign(a$proj.best[abs(a$proj.best) == max(abs(a$proj.best))])
        index <- (1:p) * (abs(a$proj.best) == max(abs(a$proj.best)))
        index <- index[index > 0]
        if (G == 2) {
            class <- i.class
        }
        else {
            m <- tapply(c(proj.data), i.class, mean)
            sd <- tapply(c(proj.data), i.class, sd)
            sd.sort <- sort.list(sd)
            m.list <- sort.list(m)
            m.sort <- sort(m)
            m.name <- as.numeric(names(m.sort))
            G <- length(m)
            dist <- 0
            split <- 0
            for (i in 1:(G - 1)) {
                if (m[m.list[i + 1]] - m[m.list[i]] > dist) {
                  split <- i
                  dist <- m[m.list[i + 1]] - m[m.list[i]]
                }
            }
            class <- rep(0, n)
            for (i in 1:split) class <- class + (i.class == m.name[i])
            class <- 2 - class
            g <- table(class)
            g.name <- as.numeric(names(g))
            G <- length(g)
            n <- nrow(i.data)
            a <- PP.optimize.anneal(PPmethod, 1, i.data, class, 
                std = TRUE, cooling, temp, energy, r, lambda)
            if (sign != sign(a$proj.best[index])) 
                a$proj.best <- -a$proj.best
            proj.data <- as.matrix(i.data) %*% a$proj.best
        }
        m.LR <- tapply(proj.data, class, mean)
        m.LR.sort <- sort(m.LR)
        LR.name <- as.numeric(names(m.LR.sort))
        var.LR <- tapply(proj.data, class, var)
        median.LR <- tapply(proj.data, class, median)
        IQR.LR <- tapply(proj.data, class, IQR)
        n.LR <- table(class)
        n.name <- as.numeric(names(n.LR))
        var.T <- sum(var.LR * n.LR)/sum(n.LR)
        if (LR.name[1] != n.name[1]) {
            temp <- n.LR[1]
            n.LR[1] <- n.LR[2]
            n.LR[2] <- temp
        }
        c1 <- (m.LR[1] + m.LR[2])/2
        c2 <- (m.LR[1] * n.LR[1] + m.LR[2] * n.LR[2])/sum(n.LR)
        max.LR <- tapply(proj.data, class, max)
        min.LR <- tapply(proj.data, class, min)
#        c3 <- sum(min.LR[2] + max.LR[1])/2
        c3 <- (m.LR[1] * var.LR[2] + m.LR[2] * var.LR[1])/sum(var.LR)        
#        c4 <- (min.LR[2] * n.LR[2] + max.LR[1] * n.LR[1])/sum(n.LR)
        c4 <- (m.LR[1] * sqrt(var.LR[2]/n.LR[2]) + m.LR[2] * sqrt(var.LR[1]/n.LR[1]))/(sqrt(var.LR[1]/n.LR[1])+sqrt(var.LR[2]/n.LR[2]))
        c5 <- (m.LR[1] * IQR.LR[2] + m.LR[2] * IQR.LR[1])/sum(IQR.LR)      
        c6 <- (m.LR[1] * (IQR.LR[2]/n.LR[2]) + m.LR[2] * (IQR.LR[1]/n.LR[1]))/((IQR.LR[1]/n.LR[1])+(IQR.LR[2]/n.LR[2]))
        C <- c(c1, c2, c3, c4,c5,c6)

        Index <- a$index.best
        Alpha <- t(a$proj.best)
        IOindexR <- NULL
        IOindexL <- NULL
        sort.LR <- as.numeric(names(sort(m.LR)))
        IOindexL <- class == sort.LR[1]
        IOindexR <- class == sort.LR[2]
        list(Index = Index, Alpha = Alpha, C = C, IOindexL = IOindexL, 
            IOindexR = IOindexR)
    }
    Tree.construct <- function(i.class, i.data, Tree.Struct, 
        id, rep, rep1, rep2, Alpha.Keep, C.Keep, PPmethod, r = NULL, 
        lambda = NULL, ...) {
        i.class <- as.integer(i.class)
        n <- nrow(i.data)
        g <- table(i.class)
        G <- length(g)
        if (length(Tree.Struct) == 0) {
            Tree.Struct <- matrix(1:(2 * G - 1), ncol = 1)
            Tree.Struct <- cbind(Tree.Struct, 0, 0, 0, 0)
        }
        if (G == 1) {
            Tree.Struct[id, 3] <- as.numeric(names(g))
            list(Tree.Struct = Tree.Struct, Alpha.Keep = Alpha.Keep, 
                C.Keep = C.Keep, rep = rep, rep1 = rep1, rep2 = rep2)
        }
        else {
            Tree.Struct[id, 2] <- rep1
            rep1 <- rep1 + 1
            Tree.Struct[id, 3] <- rep1
            rep1 <- rep1 + 1
            Tree.Struct[id, 4] <- rep2
            rep2 <- rep2 + 1
            a <- Find.proj(i.class, i.data, PPmethod, r, lambda)
            C.Keep <- rbind(C.Keep, a$C)
            Tree.Struct[id, 5] <- a$Index
            Alpha.Keep <- rbind(Alpha.Keep, a$Alpha)
            t.class <- i.class
            t.data <- i.data
            t.class <- t.class * a$IOindexL
            t.n <- length(t.class[t.class == 0])
            t.index <- sort.list(t.class)
            t.index <- sort(t.index[-(1:t.n)])
            t.class <- t.class[t.index]
            t.data <- i.data[t.index, ]
            b <- Tree.construct(t.class, t.data, Tree.Struct, 
                Tree.Struct[id, 2], rep, rep1, rep2, Alpha.Keep, 
                C.Keep, PPmethod, r, lambda)
            Tree.Struct <- b$Tree.Struct
            Alpha.Keep <- b$Alpha.Keep
            C.Keep <- b$C.Keep
            rep <- b$rep
            rep1 <- b$rep1
            rep2 <- b$rep2
            t.class <- i.class
            t.data <- i.data
            t.class <- (t.class * a$IOindexR)
            t.n <- length(t.class[t.class == 0])
            t.index <- sort.list(t.class)
            t.index <- sort(t.index[-(1:t.n)])
            t.class <- t.class[t.index]
            t.data <- i.data[t.index, ]
            n <- nrow(t.data)
            G <- length(table(t.class))
            b <- Tree.construct(t.class, t.data, Tree.Struct, 
                Tree.Struct[id, 3], rep, rep1, rep2, Alpha.Keep, 
                C.Keep, PPmethod, r, lambda)
            Tree.Struct <- b$Tree.Struct
            Alpha.Keep <- b$Alpha.Keep
            C.Keep <- b$C.Keep
            rep <- b$rep
            rep1 <- b$rep1
            rep2 <- b$rep2
        }
        list(Tree.Struct = Tree.Struct, Alpha.Keep = Alpha.Keep, 
            C.Keep = C.Keep, rep = rep, rep1 = rep1, rep2 = rep2)
    }
    C.Keep <- NULL
    Alpha.Keep <- NULL
    Tree.Struct <- NULL
    id <- 1
    rep1 <- 2
    rep2 <- 1
    rep <- 1
    if (PPmethod == "LDA" && weight) 
        method <- 1
    else if (PPmethod == "LDA" && !weight) 
        method <- 2
    else if (PPmethod == "Lp") 
        method <- 3
    else if (PPmethod == "Gini") 
        method <- 4
    else if (PPmethod == "Ent") 
        method <- 5
    else if (PPmethod == "PDA") 
        method <- 6
    else stop("Wrong PPmethod")
    Tree.final <- Tree.construct(i.class, i.data, Tree.Struct, 
        id, rep, rep1, rep2, Alpha.Keep, C.Keep, PPmethod, r, 
        lambda)
    Tree.Struct <- Tree.final$Tree.Struct
    colnames(Tree.Struct)<-c("id","L.node.ID","R.F.node.ID","Coef.ID","Index")
    Alpha.Keep <- Tree.final$Alpha.Keep
    C.Keep <- Tree.final$C.Keep
    list(Tree.Struct = Tree.Struct, Alpha.Keep = Alpha.Keep, 
        C.Keep = C.Keep)
}


PP.classify <- function(test.data, true.class, Tree.result, Rule, ...) {
    test.data<-as.matrix(test.data)
    true.class<-as.matrix(true.class); 
    if(nrow(true.class)==1) true.class<-t(true.class)
    if (!is.numeric(true.class)) {
        class.name<-names(table(true.class))
        temp<-rep(0,nrow(true.class))
        for(i in 1:length(class.name))
            temp<-temp+(true.class==class.name[i])*i
        true.class<-temp
    } 
#############
    PP.Classification <- function(Tree.Struct, test.class.index, IOindex,
                                  test.class, id, rep) {
        if (Tree.Struct[id,4] == 0) {
            i.class <- test.class
            i.class[i.class > 0] <- 1
            i.class <- 1 - i.class
            test.class <- test.class + IOindex * i.class * Tree.Struct[id, 3]
            return(list(test.class=test.class, rep=rep))
        } else {  
            IOindexL <- IOindex * test.class.index[rep,]
            IOindexR <- IOindex * (1 - test.class.index[rep,])
            rep <- rep + 1
            a <- PP.Classification(Tree.Struct, test.class.index, IOindexL,
                                   test.class, Tree.Struct[id,2], rep)
            test.class <- a$test.class
            rep <- a$rep;
            a <- PP.Classification(Tree.Struct, test.class.index, IOindexR,
                                   test.class, Tree.Struct[id,3], rep)
            test.class <- a$test.class
            rep <- a$rep
        }
        list(test.class=test.class, rep=rep)
    }
##############
    
    PP.Class.index <- function(class.temp, test.class.index, test.data,
                               Tree.Struct, Alpha.Keep, C.Keep, id,Rule) {
        class.temp <- as.integer(class.temp)
        if (Tree.Struct[id,2] == 0) {
            return(list(test.class.index=test.class.index,
                        class.temp=class.temp))
        } else {
            t.class <- class.temp 
            t.n <- length(t.class[t.class == 0])
            t.index <- sort.list(t.class)
            if (t.n) t.index <- sort(t.index[-(1:t.n)])
            t.data <- test.data[t.index,]
            id.proj <- Tree.Struct[id,4]
            
            proj.test <- as.matrix(test.data) %*%
                as.matrix(Alpha.Keep[id.proj,])
            ##  proj.test<-(proj.test-mean(proj.test))
            proj.test <- as.double(proj.test)
            class.temp <- t(proj.test < C.Keep[id.proj,Rule]) 
            test.class.index <- rbind(test.class.index, class.temp)
            a <- PP.Class.index(class.temp, test.class.index, test.data,
                                Tree.Struct, Alpha.Keep, C.Keep,
                                Tree.Struct[id,2], Rule)
            test.class.index <- a$test.class.index
            a<-PP.Class.index(1 - class.temp, test.class.index, test.data,
                              Tree.Struct, Alpha.Keep, C.Keep,
                              Tree.Struct[id,3], Rule)
            test.class.index <- a$test.class.index;
        }
        list(test.class.index=test.class.index, class.temp=class.temp)
    }
    
    
#############
    n <- nrow(test.data)
    class.temp <- rep(1, n)
    test.class.index <- NULL
    temp <- PP.Class.index(class.temp, test.class.index, test.data,
                           Tree.result$Tree.Struct, Tree.result$Alpha.Keep,
                           Tree.result$C.Keep, 1, Rule)
    test.class <- rep(0, n)
    IOindex <- rep(1, n)
    rep <- 1
    temp <- PP.Classification(Tree.result$Tree.Struct, temp$test.class.index,
                              IOindex, test.class, 1, 1)
    predict.error <- sum(true.class != temp$test.class)
    predict.class <- temp$test.class
    list(predict.error=predict.error, predict.class=predict.class)
}

