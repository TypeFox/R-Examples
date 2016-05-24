
#############################################
#
# PP index calculation
#
#############################################
PPindex.class <- function(PPmethod, data, class, weight=TRUE, r=NULL,
                          lambda=NULL, ...) {
    if (PPmethod =="LDA") 
        index <- PPindex.LDA(data, class, weight)
    else if (PPmethod =="Lp") 
        index <- PPindex.Lp(data,class,r)
    else if (PPmethod =="PDA") 
        index <- PPindex.PDA(data,class,lambda)
    else
        stop("You need to select PPmethod!")
    index
}

PPindex.LDA <- function(data, class, weight=TRUE, ...) {
    data <- as.matrix(data)
    class <- as.matrix(class)
    if (ncol(class) != 1) class <- t(class)
    p <- ncol(data)
    n <- nrow(data)
    ngroup <- table(class)
    groups <- length(ngroup)
    
    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    class.i <- matrix(class.num[as.character(class)], ncol=1)

    val <- 1
    class.sort <- sort.list(class.i)
    class <- class.i[class.sort]
    data <- data[class.sort, ]
    ngroup <- table(class.i)
    groups <- length(ngroup)
    gname <- names(ngroup)
  
    LDA <- if (weight)
        .C("discriminant1", 
           as.integer(n),
           as.integer(p),
           as.integer(groups),
           x=as.double(data),
           class=as.integer(class), 
           as.integer(gname),
           ngroup = as.integer(ngroup),
           index = as.double(val),
           PACKAGE="PPtree")
    else
        .C("discriminant2", 
           as.integer(n),
           as.integer(p),
           as.integer(groups),
           x=as.double(data),
           class=as.integer(class), 
           as.integer(gname),
           ngroup = as.integer(ngroup),
           index = as.double(val),
           PACKAGE="PPtree")
    LDA$index
}

PPindex.PDA <- function(data, class, lambda, ...) {
    if (is.null(lambda)) 
        stop("You need to use parameter lambda!")
    data <- as.matrix(data)
    class <- as.matrix(class)
    if (ncol(class) != 1) class <- t(class)
    p <- ncol(data)
    n <- nrow(data)
    ngroup <- table(class)
    groups <- length(ngroup)

    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    class.i <- matrix(class.num[as.character(class)], ncol=1)

    val<-1
    class.sort <- sort.list(class.i)
    class <- class.i[class.sort]
    data <- data[class.sort,]
    ngroup <- table(class.i)
    groups <- length(ngroup)
    gname <- names(ngroup)
    
    PDA <- .C("pda", 
              as.integer(n),
              as.integer(p),
              as.integer(groups),
              x=as.double(data),
              class=as.integer(class), 
              as.integer(gname),
              ngroup = as.integer(ngroup),
              index = as.double(val),
              as.double(lambda),
              PACKAGE="PPtree")
    PDA$index
}

PPindex.Lp <- function(data, class, r, ...) {
    if (is.null(r)) 
        stop("You need to select parameter r")
    
    data<-as.matrix(data);
    class<-as.matrix(class);
    if(ncol(class)!=1) class<-t(class)
    p<-ncol(data);n<-nrow(data);
    ngroup<-table(class)
    groups<-length(ngroup)
    
    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    class.i <- matrix(class.num[as.character(class)], ncol=1)
    
    val <- 1
    class.sort <- sort.list(class.i)
    class <- class.i[class.sort]
    data <- data[class.sort,]
    ngroup <- table(class.i)
    groups <- length(ngroup)
    gname <- names(ngroup)
    
    LDA <- .C("discriminant3", 
              as.integer(n),
              as.integer(p),
              as.integer(groups),
              x=as.double(data),
              class=as.integer(class), 
              as.integer(gname),
              ngroup=as.integer(ngroup),
              index=as.double(val),
              as.integer(r),
              PACKAGE="PPtree")
    LDA$index
}

##########################################################
PP.optimize.random <- function(PPmethod, projdim, data, class, std=TRUE,
                               cooling=0.99, temp=1, r=NULL, lambda=NULL,
                               weight=TRUE, ...) {
    data <- as.matrix(data)
    class <- as.matrix(class)
    pp <- ncol(data)
    if (std) {
        remove <- (1:pp)*(apply(data,2,sd) == 0)
        remove <- remove[remove != 0]
        if (length(remove) != 0) {
            data <- data[, -remove]
        } 
        data <- scale(data) # apply(data,2,function(x){(x-mean(x))/sd(x)})
    }

    if (ncol(class) != 1) class <- t(class)
    p <- ncol(data)
    n <- nrow(data)
    ngroup <- table(class)
    groups <- length(ngroup)

    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    class.i <- matrix(class.num[as.character(class)], ncol=1)
    val <- 1
    proj <- matrix(1, p, projdim)
    class.sort <- sort.list(class.i)
    class <- class.i[class.sort]
    data <- data[class.sort, ]
    ngroup <- table(class.i)
    groups <- length(ngroup)
    gname <- names(ngroup)

    if (PPmethod =="LDA" && weight==TRUE) 
        method <- 1
    else if(PPmethod =="LDA" && weight==FALSE) 
        method <- 2
    else if (PPmethod =="Lp") {
        method <- 3
        if (is.null(r)) 
            stop("You need to select parameter r!")
    }
    else if (PPmethod =="Gini") 
        method <- 4
    else if (PPmethod =="Ent") 
        method <- 5
    else  if(PPmethod == "PDA") {
        method <- 6
        if (is.null(lambda)) 
            stop("You need to select parameter lambda !")
    }
    else stop("You need to select PPmethod!")
    
    temp.start <- 1
    temp.end <- 0.001
    heating <- 1
    restart <- 1
    success <- 0
    maxproj <- 1000  
    
    Opt <- .C("optimize1", 
              as.integer(n),
              as.integer(p),
              as.integer(groups),
              x=as.double(data),
              as.integer(class),
              as.integer(gname), 
              as.integer(ngroup),
              as.integer(method),
              cooling=as.double(cooling),
              temp=as.double(temp),
              projdim=as.integer(projdim),
              index=as.double(val),
              proj=as.double(proj),
              as.integer(r),
              as.double(lambda),
              PACKAGE="PPtree")
    index.best <- Opt$index
    if (pp != p) {
        proj.best <- matrix(0, pp, projdim)
        proj.best[-remove, ] <- matrix(Opt$proj, ncol=projdim)
    } else {
        proj.best <- matrix(Opt$proj, ncol=projdim)
    }
    list(index.best=index.best, proj.best=proj.best)
}

#########################
PP.optimize.anneal <- function(PPmethod, projdim, data, class, std=TRUE,
                               cooling=0.999, temp=1, energy=0.01, r=NULL,
                               lambda=NULL, weight=TRUE, ...) {
    data <- as.matrix(data)
    class <- as.matrix(class)
    pp <- ncol(data)
    if (std) {
        remove <- (1:pp)*(apply(data,2,sd) == 0)
        remove <- remove[remove != 0]
        if (length(remove)!=0) {
            warning("Removing variables with 0 variance")
            data<-data[,-remove]
        }
        data <- scale(data)
    }

    if (ncol(class) != 1) class <- t(class)
    p <- ncol(data)
    n <- nrow(data);
    ngroup <- table(class)
    groups <- length(ngroup)
    
    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    class.i <- matrix(class.num[as.character(class)], ncol = 1)

    val <- 1
    proj <- matrix(1, p, projdim)
    class.sort <- sort.list(class.i)
    class <- class.i[class.sort]
    data <- data[class.sort,]
    ngroup <- table(class.i)
    groups <- length(ngroup)
    gname <- names(ngroup)
    if(PPmethod =="LDA" && weight) 
        method <- 1
    else if (PPmethod =="LDA" && !weight) 
        method <- 2
    else if (PPmethod == "Lp") {
        method <- 3
        if (is.null(r)) 
            stop("You need to select parameter r!")
    }
    else if (PPmethod == "Gini") 
        method <- 4
    else if(PPmethod == "Ent") 
        method <- 5
    else if (PPmethod == "PDA") {
        method<- 6
        if (is.null(lambda)) 
            stop("You need to select parameter lambda !")
    }
    else stop("You need to select PPmethod!")

    Opt <- .C("optimize2", 
              as.integer(n),
              as.integer(p),
              as.integer(groups),
              x=as.double(data),
              as.integer(class), 
              as.integer(gname),
              as.integer(ngroup),
              as.integer(method),
              cooling=as.double(cooling),
              temp=as.double(temp),
              projdim=as.integer(projdim),
              index=as.double(val),
              proj=as.double(proj),
              as.double(energy),
              as.integer(r),
              as.double(lambda),
              PACKAGE="PPtree")
    index.best<-Opt$index
    if (pp != p) {
        proj.best <- matrix(0, pp, projdim)
        proj.best[-remove,] <- matrix(Opt$proj, ncol=projdim)
    } else {
        proj.best <- matrix(Opt$proj, ncol=projdim)
    }
    list(index.best=index.best, proj.best=proj.best)
}

###############
PP.optimize.Huber <- function(PPmethod, projdim, data, class, std=TRUE,
                              cooling=0.99, temp=1, r=NULL, lambda=NULL,
                              weight=TRUE, ...) {
    data <- as.matrix(data)
    class <- as.matrix(class);
    pp <- ncol(data)
    if (std) {
        remove <- (1:pp)*(apply(data,2,sd) == 0)
        remove <- remove[remove != 0]
        if (length(remove) != 0) {
            warning("Remove variables with 0 variance")
            data<-data[, -remove]
        }
        data <- scale(data)
    }
    
    if (ncol(class) != 1) class <- t(class)
    p <- ncol(data)
    n <- nrow(data);
    ngroup <- table(class)
    groups <- length(ngroup)
    
    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    class.i <- matrix(class.num[as.character(class)], ncol=1)
    
    val <- 1
    proj <- matrix(1, p, projdim)
    class.sort <- sort.list(class.i)
    class <- class.i[class.sort]
    data <- data[class.sort,]
    ngroup <- table(class.i)
    groups <- length(ngroup)
    gname <- names(ngroup)
    if (PPmethod =="LDA" && weight==TRUE) 
        method <- 1
    else if (PPmethod == "LDA" && weight==FALSE) 
        method <- 2
    else if (PPmethod == "Lp") {
        method <- 3
        if (is.null(r)) stop("You need to select parameter r !")
    }
    else if (PPmethod == "Gini") 
        method <- 4
    else if (PPmethod == "Ent") 
        method <- 5
    else if (PPmethod == "PDA") {
        method <- 6
        if (is.null(lambda)) stop("You need to select parameter lambda !")
    }
    else stop("You need to select PPmethod!")

    Opt <- .C("optimize3", 
              as.integer(n),
              as.integer(p),
              as.integer(groups),
              x=as.double(data),
              as.integer(class), 
              as.integer(gname),
              as.integer(ngroup),
              as.integer(method),
              cooling=as.double(cooling),
              temp=as.double(temp),
              projdim=as.integer(projdim),
              index=as.double(val),
              proj=as.double(proj),
              as.integer(r),
              PACKAGE="PPtree")
    index.best <- Opt$index
    if (pp !=p) {
        proj.best <- matrix(0, pp, projdim)
        proj.best[-remove,] <- matrix(Opt$proj, ncol=projdim)
    } else {
        proj.best <- matrix(Opt$proj, ncol=projdim)
    }
    list(index.best=index.best, proj.best=proj.best)
}

#########################
PP.optimize.plot <- function(PP.opt, data, class, std=TRUE) { 
    pp<-ncol(data)
    n<-nrow(data);
    if(std) {
        remove <- c(1:pp)*(apply(data,2,sd) == 0)
        remove <- remove[remove != 0]
        if (length(remove) != 0) {
            warning("Remove variables with 0 variance")
            data <- data[,-remove]
        } 
        data.s <- scale(data)
        if (pp != ncol(data.s)) {
            data <- matrix(0, n, pp)
            data[,-remove] <- matrix(data.s, ncol=pp)
        }
    }
    data<-as.matrix(data)
    n<-nrow(data)
    class <- as.matrix(class)
    if (ncol(class) != 1) class <- t(class)
    ngroup <- table(class)
    groups <- length(ngroup)

    class.num <- 1:groups
    names(class.num) <- names(ngroup)
    t.class <- matrix(class.num[as.character(class)], ncol=1)

    p <- ncol(PP.opt$proj.best)
    proj.data <- as.matrix(data) %*% PP.opt$proj.best
    labels <- NULL
    for (i in 1:p)
        labels <- c(labels, paste("PP",i,sep=""))
    if (p == 1) { 
        hist(proj.data, main=" ", xlab=labels[1]) 
        ngroup <- table(t.class)
        ng <- length(ngroup)
        sort.index <- sort.list(t.class)
        sort.data <- sort(t.class)
        k <- 0
        step <- n / 10 / ng
        for (j in 1:ng) {
            ngg <- ngroup[j]
            for(i in 1:ngg) {
                k <- k + 1
                text(proj.data[sort.index[k]], step*j, 
                     as.character(t.class[sort.index[k]]), cex=2)  
            }
        }
    } else {
        if (p == 2) {
            par(pty='s')
            plot(proj.data[,1], proj.data[,2], type='n', xlab=labels[1],
                 ylab=labels[2])
            text(proj.data[,1], proj.data[,2], as.character(t.class), cex=2) 
        } else {
            g <- length(table(t.class))
            t.class <- as.numeric(t.class)
            color <- colors()[(1:g) * 10]
            pairs(proj.data, pch=21, bg=color[t.class], labels)
        }
    }
}