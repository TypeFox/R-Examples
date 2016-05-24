## ===============================================
## Plotting individual trajectories of a categorical
## longitudinal variable
## ===============================================

ctplot <- function(x,y,id,weights=NULL,
                   cov=NULL,subset=NULL,type="distinctive",
                   embedding="most-frequent",x.order=FALSE,
                   x.orderalign=ifelse(split=="last","last","first"),
                   y.optimalphabet=FALSE,R=NULL,
                   main="",sub=NULL,mtext=TRUE,
                   xlab="order position",ylab="",xlim=NULL,ylim=NULL,
                   grid.col="white",grid.fill="grey90",
                   grid.scale=1/4,grid.shape="default",grid.lwd=0,
                   cpal=NULL,alpha=1,lcourse="upwards",lorder="background",
                   lweights=TRUE,lwd.min=0.5,lwd.max=2,lty=1,pch=4,
                   cex=1,border=grid.col,border.lwd=grid.lwd/2,
                   sf.cex=1,sf.cex.leaves=1,sf.pch=16,rotate=FALSE,
                   split=NULL,layout=NULL,
                   show.type=1,show=c(0,1),hide.col="grey75",
                   hide.s=0.1,hide.na.cost=1,
                   print=FALSE,return.data=FALSE,...)
  {

    if (!is.factor(id)) {
      id <- factor(id)
    }

    ## ======
    ## checks
    ## ======

    if (!is.null(weights)) {
      if (length(weights)!=nlevels(id)) {
        stop("[!] impossible to link id and weights vector. See help")
      } else {
        names(weights) <- levels(id)
      }
    }

    if (!is.null(cov)) {
      if (length(cov)!=nlevels(id)) {
        stop("[!] impossible to link id and cov vector. See help")
      } else {
        names(cov) <- levels(id)
      }
    }

    if (!x.order&!is.null(split)) {
      warning(" [!] split option is reserved for ordered data. ignore assigned split argument")
      split <- NULL
    }

    if (sum(c(length(x),length(y))!=length(id))!=0) {
      stop("[!] input vectors x, y and id have different lengths")
    }

    if (!(type%in%c("jitter","distinctive","non-embeddable"))) {
      stop("[!] invalid trajectory type")
    }

    if (!(embedding%in%c("most-frequent","uniformly"))) {
      stop("[!] invalid embedding type")
    }

    if (!lweights) {
      lwd.max <- lwd.min
    }

    if ((min(show)<0)|(max(show)>1)) {
      stop("[!] invalid show parameters")
    }

    if (is.numeric(y)) {
      y <- cut(y,5,ordered_result=TRUE)
      warning("[!] transformed y to arbirary ordinal scale")
    }

    if (is.ordered(y)) {
      ytype <- "ordinal"
    } else {
      ytype <- "categorical"
    }

    ## ==================
    ## preparing the data
    ## ==================

    if (print) cat(" [>] preparing the data\n")

    ## extract subset
    ## ==============

    if (!is.null(subset)) {
      if (length(subset)>length(unique(subset))) {
        warning("[!] duplicated ids in subset vector")
        subset <- unique(subset)
      }
      subset <- as.character(subset)
      subset.LONG <- id%in%subset
      x <- x[subset.LONG]
      y <- y[subset.LONG]
      id <- id[subset.LONG]
      id <- factor(id,levels=subset)
      if (!is.null(weights)) {
        subset.ID <- names(weights)%in%subset
        weights <- weights[subset.ID]
      }
      if (!is.null(cov)) {
        subset.ID <- names(cov)%in%subset
        cov <- cov[subset.ID]
      }
    }

    ## erase useless NA entries
    ## ========================

    ## NAs in weights
    if (!is.null(weights)) {
      subscripts.ID <- is.na(weights)
      subscripts.LONG <- id%in%names(weights[subscripts.ID])
      if (sum(subscripts.LONG)>0) {
        subscripts.id <- setdiff(names(weights),
                                 names(weights[!subscripts.ID]))
        warning(paste(" [!] ignore NA entries of ids: ",
                      paste(subscripts.id,collapse=", "),
                      " in weights vector",sep=""))
        weights <- weights[!subscripts.ID]
        id <- factor(id[!subscripts.LONG],levels=names(weights))
        x <- x[!subscripts.LONG]
        y <- y[!subscripts.LONG]
        if (!is.null(cov)) {
          cov <- cov[!subscripts.ID]
        }
      }
    }

    ## NAs in covariate
    if (!is.null(cov)) {
      subscripts.ID <- is.na(cov)
      subscripts.LONG <- id%in%names(cov[subscripts.ID])
      if (sum(subscripts.LONG)>0) {
        subscripts.id <- setdiff(names(cov),names(cov[!subscripts.ID]))
        warning(paste(" [!] ignore NA entries of ids: ",
                      paste(subscripts.id,collapse=", "),"in cov",sep=""))
        cov <- cov[!subscripts.ID]
        weights <- weights[!subscripts.ID]
        id <- factor(id[!subscripts.LONG],levels=names(cov))
        x <- x[!subscripts.LONG]
        y <- y[!subscripts.LONG]
      }
    }

    ## NAs in x, y or id
    subscripts.LONG <- is.na(id)|is.na(x)|is.na(y)
    subscripts.id <- intersect(id,id[subscripts.LONG])
    if (sum(subscripts.LONG)>0) {
      warning(paste(" [!] ignore NA entries of",
                    paste(subscripts.id,collapse=", ")))
      id <- id[!subscripts.LONG]
      x <- x[!subscripts.LONG]
      y <- y[!subscripts.LONG]
    }

    ## standardize weights
    if (!is.null(weights)) {
      if (!is.null(cov)) {
        weights.cov <- table(cov)[as.integer(cov)]
        weights.sum <- tapply(X=weights,INDEX=cov,FUN=sum)[as.integer(cov)]
        weights <- weights/weights.sum*weights.cov
        names(weights) <- levels(id)
      } else {
        weights <- weights/sum(weights)
      }
    }

    ## check the x variable
    ## ====================

    if (!x.order) {
      if (is.numeric(x)) {
        if (length(unique(x))>1) {
          xdiff <- diff(sort(unique(x)))
          if (sum((xdiff/min(xdiff)-round(xdiff/min(xdiff),0))!=0)==0) {
            x <- factor(x,levels=seq(min(x),max(x),min(xdiff)))
          } else {
            warning("[!] Problems with distances between order positions. The order positions will not be illustrated adequately.")
            x <- factor(x)
          }
        } else {x <- factor(x)}
      }
      xtickat <- 1:nlevels(x)
      xticklab <- levels(x)
      x <- as.integer(x)
    }

    ## weight variable
    ## ===============

    if (is.null(weights)) {
      weights <- rep(1,nlevels(id))
      names(weights) <- levels(id)
    }

    ## check the id variable
    ## =====================

    ## ids without y observations
    usedid <- intersect(levels(id),unique(id))
    nid.tot <- nlevels(id)
    id <- factor(id,levels=usedid)
    weights.use <- weights[names(weights)%in%usedid]
    nid.use <- nlevels(id)
    id.int <- as.integer(id)

    ## check the y variable
    ## ====================

    ny <- nlevels(y)
    y.int <- as.integer(y)

    ## order the data entries
    ## ======================

    ord <- order(id,x,y)
    id <- id[ord]
    id.int <- id.int[ord]
    x <- x[ord]
    y <- y[ord]
    y.int <- y.int[ord]

    ## check covariate
    ## ===============

    if (!is.null(cov)) {
      cov <- factor(cov)
      nid.cov.tot <- table(cov)
      cov.old <- cov
      cov <- cov[names(cov)%in%usedid]
      nid.cov.use <- table(cov)
    } else {
      ncov <- 1
      cov.int <- rep(1,nid.use)
      nid.cov.tot <- nid.tot
      nid.cov.use <- nid.use
    }
    if (!is.null(cov)) {
      if (!is.factor(cov)) {
        cov <- factor(cov)
      }
      if (sum(nid.cov.tot==0)>0) {
        warning(paste(" [!] the plot will omit the following covariate categories:",paste(names(nid.cov.tot)[nid.cov.tot==0],collapse=", "),sep=" "))
        cov <- factor(cov,levels=names(nid.cov.tot)[nid.cov.tot>0])
      }
      ncov <- nlevels(cov)
      cov.int <- as.integer(cov)
    } else {
      ncov <- 1
      cov.int <- rep(1,nid.use)
    }

    ## build orders
    ## ============

    if (x.order) {
      x <- unlist(tapply(X=x,INDEX=list(id),FUN=retorders))
      xtickat <- xticklab <- 1:max(x)
    }

    ## consider order alignment
    ## ========================

    if (x.order) {
      if (is.character(x.orderalign)) {
        if (x.orderalign=="first") {
          xtickat <- xticklab <- 1:max(x)
        } else {
          if (x.orderalign=="last") {
            align <- tapply(X=x,INDEX=id.int,FUN=max)
            x <- x-align[id.int]
            x <- x-min(x)+1
            xtickat <- 1:max(x)
            xticklab <- seq(1-max(x),0)
          } else {
            if (x.orderalign%in%levels(y)) {
              subscripts <- y==x.orderalign
              subscripts <- unlist(tapply(X=subscripts,
                                          INDEX=list(id.int),
                                          FUN=remfolltrues))
              aglevs <- c("no align",paste(x.orderalign,"align"))
              align <- rep(1,nid.use)
              tab <- table(y,id)[x.orderalign,]
              align[which(tab>0)] <- x[subscripts]
              aligngroup <- factor(aglevs[(tab>0)+1])
              if (!is.null(cov)) {
                cov <- factor(paste(cov,aligngroup,sep="/"))
              } else {
                cov <- aligngroup
              }
              ncov <- nlevels(cov)
              cov.int <- as.integer(cov)
              x <- x-align[id.int]
              x <- x-min(x)+1
              xtickat <- 1:max(x)
              xticklab <- xtickat-x[which(y==x.orderalign)[1]]
            }
          }
        }
      }
    }
    nx <- max(x)

    ## ## construct x and y lists (one element per id)
    ## ## ============================================

    y.list <- tapply(X=y.int,
                    INDEX=list(id.int),
                    FUN=function(x){return(x)})
    x.list <- tapply(X=x,
                    INDEX=list(id.int),
                    FUN=function(x){return(x)})

    ## construct trajectory strings (one element per id)
    ## =================================================

    traj.seqe <- seqecreate(id=id.int,timestamp=x,event=y.int)
    traj.string <- tapply(X=y.int,INDEX=list(id.int,x),paste,collapse=",")
    traj.string <- apply(X=traj.string,MARGIN=1,FUN=createstring)
    traj.order <- order(nchar(traj.string),decreasing=TRUE)

    ## ===========================================
    ## find distinctive trajectories
    ## ===========================================

    if (print) cat(" [>] find distinctive trajectories\n")

    traj.group <- c()
    traj.id <- c()
    group.id <- c()
    group.name <- c()

    embedded <- c()
    if (type=="non-embeddable") {
      while (sum(!(1:nid.use)%in%embedded)>0) {
        i <- which(!traj.order%in%embedded)[1]
        gn <- ifelse(length(traj.group)==0,1,max(traj.group)+1)
        group.id <- append(group.id,traj.order[i])
        group.name <- append(group.name,gn)

        ## check for embeddable trajectories
        if (is.null(embedded)) {
          subscripts <-
            which(unlist(lapply(X=traj.string,
                                FUN=grepl,traj.string[traj.order[i]],
                                fixed=TRUE)))
        } else {
          subscripts <-
            which(unlist(lapply(X=traj.string[-traj.order[1:(i-1)]],
                                FUN=grepl,traj.string[traj.order[i]],
                                fixed=TRUE)))
          subscripts <- (1:nid.use)[-traj.order[1:(i-1)]][subscripts]
        }
        traj.id <- append(traj.id,subscripts)
        traj.group <- append(traj.group,rep(gn,length(subscripts)))
        embedded <- unique(c(embedded,subscripts))
      }
    } else if (type=="distinctive") {
      traj.group <- as.integer(factor(traj.string))
      traj.id <- 1:nid.use
      group.name <- 1:max(traj.group)
      group.id <- tapply(traj.id,traj.group,function(x){x[1]})
    } else if (type=="jitter") {
      traj.group <- 1:length(traj.string)
      traj.id <- 1:nid.use
      group.name <- 1:max(traj.group)
      group.id <- 1:nid.use
    }

    ## erase difficult embeddings difficult to interprete
    ## in case a covariate vector was entered
    ## ==================================================

    tab <- table(traj.id,traj.group)
    if ((type=="non-embeddable")&(ncov>1)) {
      uemb <- rowSums(tab)==1
      for (i in 1:ncov) {
        ## for which non-embeddedable trajectories there are observed
        ## trajectories that were solely embedded into it

        uemb.groups <- traj.group[traj.id%in%which(uemb&(cov.int==i))]
        for (j in which((cov.int==i)&(!uemb)))
          ## erase embeddings
          if (sum(uemb.groups%in%traj.group[traj.id%in%j])>0) {
            subscripts <- (traj.id%in%j)&(!traj.group%in%uemb.groups)
            traj.id <- traj.id[!subscripts]
            traj.group <- traj.group[!subscripts]
          }
      }
    }

    ## construct the plot data
    ## =======================

    pointdata <- data.frame(id=id.int[id.int%in%group.id],
                            x=x[id.int%in%group.id],
                            y=y.int[id.int%in%group.id])
    pointdata$group <- as.integer(factor(pointdata$id,
                                         levels=group.id,
                                         labels=group.name))

    ## if lines are demanded to run from the highest to lowest y category
    ## in case of simultaneous observations
    if (lcourse=="downwards") {
      pointdata <-
        pointdata[order(pointdata$group,pointdata$x,
                        factor(pointdata$y,
                               levels=rev(sort(unique(pointdata$y))))),]
    }

    ## find multiple equal simultaneous observations
    ## and note their frequency
    finddoubles <- data.frame(table(x=pointdata$x,
                                    y=pointdata$y,
                                    group=pointdata$group))
    finddoubles <- finddoubles[finddoubles$Freq!=0,]
    pointdata <- merge(x=unique(pointdata),
                       y=finddoubles,
                       by=c("x","y","group"),
                       sort=FALSE)

    ## extract coordinate pairs for plotting the lines using
    ## segments (necessary to allow different line widths)
    ## =====================================================

    ngroup <- max(pointdata$group)
    linedata <- NULL

    for (g in 1:ngroup) {
      subs <- which(pointdata$group==g)
      nsubs <- length(subs)
      if (length(subs)>1) {
        linedata.tmp <-
          data.frame(group=rep(g,nsubs-1),
                     x0=pointdata$x[subs[-nsubs]],
                     y0=pointdata$y[subs[-nsubs]],
                     x1=pointdata$x[subs[-1]],
                     y1=pointdata$y[subs[-1]])
        linedata <- rbind(linedata,linedata.tmp)
      }
    }

    ## ===============================================
    ## calculate weighted trajectory shape group sizes
    ## ===============================================

    pointdata[,paste("n",1:ncov,sep=".")] <- 1
    if (!is.null(linedata)) {
      linedata[,paste("n",1:ncov,sep=".")] <- 0
      linedata[,paste("lwd",1:ncov,sep=".")] <- 1
    }
    if (print) cat(" [>] determine symbol sizes and line widths\n")

    tab <- table(traj.id,traj.group)
    w <- 1/rowSums(tab)
    minx <- unlist(lapply(X=x.list,FUN=min))
    maxx <- unlist(lapply(X=x.list,FUN=max))

    ## option: trajectories that are embeddable into multiple
    ## non-embeddable trajectories are embedded into the most
    ## frequent non-embeddable trajectory

    if (embedding=="most-frequent") {

      weights.group <- tapply(X=weights.use[w==1],
                              INDEX=list(traj.group[traj.id%in%which(w==1)],
                                factor(cov.int[w==1],levels=1:ncov)),
                              FUN=sum)
      weights.group[is.na(weights.group)] <- 0

      ## perform unique group assignements
      for (i in which(w!=1)) {
        group.iin <- traj.group[traj.id==i]
        group.iin <-
          group.iin[which.max(weights.group[group.iin,cov.int[i]])]
        subs <- which((traj.group!=group.iin)&(traj.id==i))
        traj.id <- traj.id[-subs]
        traj.group <- traj.group[-subs]
      }
      w <- rep(1,length(w)) # set all
    }

    ## point weights
    for (i in 1:nrow(pointdata)) {
      if (i==1) {
        isubs <- traj.id[traj.group==pointdata$group[i]]
      } else {
        if (pointdata$group[i]!=pointdata$group[i-1]) {
          isubs <- traj.id[traj.group==pointdata$group[i]]
        }
      }
      pointdata[i,paste("n",1:ncov,sep=".")] <-
        unlist(tapply(X=weights.use[isubs]*w[isubs]*
                      ((minx[isubs]<=pointdata$x[i])&
                       (maxx[isubs]>=pointdata$x[i])),
                      INDEX=list(factor(cov.int[isubs],levels=1:ncov)),
                      FUN=sum))
    }
    pointdata[is.na(pointdata)] <- 0

    ## line weights
    if (!is.null(linedata)) {
      for (i in 1:nrow(linedata)) { # loop over all linedata entries
        if (i==1) {
          isubs <- traj.id[traj.group==linedata$group[i]]
        } else {
          if (linedata$group[i]!=linedata$group[i-1]) {
            isubs <- traj.id[traj.group==linedata$group[i]]
          }
        }
        linedata[i,paste("n",1:ncov,sep=".")] <-
          unlist(tapply(X=weights.use[isubs]*w[isubs]*
                        ((minx[isubs]<=linedata$x0[i])&
                         (maxx[isubs]>=linedata$x1[i])),
                        INDEX=list(factor(cov.int[isubs],levels=1:ncov)),
                        FUN=sum))
      }
      linedata[is.na(linedata)] <- 0
    }

    ## line widths
    if (!is.null(linedata)) {
      if (!is.null(cov)) {
        maxprop <-
          apply(as.matrix(linedata[,paste("n",1:ncov,sep=".")],
                          ncol=ncov),2,max,na.rm=TRUE)/nid.cov.tot

        facm <- diag(maxprop/max(maxprop))
        scalem <- apply(as.matrix(linedata[,paste("n",1:ncov,sep=".")],
                                  ncol=ncov),2,max,na.rm=TRUE)
        scalem[scalem==0] <- 1
        scalem <- scale(x=linedata[,paste("n",1:ncov,sep=".")],
                        center=FALSE,scale=scalem)
        scalem <- scalem%*%facm
        linedata[,paste("lwd",1:ncov,sep=".")] <- as.data.frame(scalem)
      } else {
        linedata[,paste("lwd",1:ncov,sep=".")] <-
          linedata[,paste("n",1:ncov,sep=".")]/
            max(linedata[,paste("n",1:ncov,sep=".")])
      }
    }

    ## ========================
    ## optimize order of y axis
    ## ========================

    if (y.optimalphabet)
      {
        if (print) cat(" [>] optimize order\n")

        permn <- permn(1:nlevels(y))
        if (is.numeric(R)) {
          permn <- permn[sample(1:length(permn),R)]
        }

        dist <- rep(0,length(permn))
        for (i in 1:length(permn))
          {
            pointdata$y.tmp <- as.integer(factor(x=pointdata$y,
                                                 levels=permn[[i]]))
            for (g in 1:ngroup)
              {
                IND <- pointdata$group==g
                dist[i] <- dist[i]+sqrt(sum(diff(pointdata$y.tmp[IND])^2+
                                            diff(pointdata$x[IND])))
              }
          }
        pointdata$y <- as.integer(factor(x=pointdata$y,
                                         levels=permn[[which.min(dist)]]))
        y <- factor(y,levels=levels(y)[permn[[which.min(dist)]]])
      }

    ## ==============================
    ## Find plot positions for groups
    ## ==============================

    if (print) cat(" [>] determine plot coordinates\n")

    maxgs <- apply(X=as.matrix(x=pointdata[,paste("n",1:ncov,sep=".")],
                     ncol=ncov),MARGIN=2,FUN=tapply,
                   list(pointdata$group),max)
    maxgs <- scale(x=matrix(maxgs,ncol=ncov),
                   center=FALSE,scale=nid.cov.tot)
    maxgs.tot <- apply(X=as.matrix(x=maxgs,ncol=ncov),MARGIN=1,FUN=max)
    ord <- order(maxgs.tot,decreasing=TRUE)

    ## arrangement with non-overlapping positions
    ## ==========================================

    ## create initial grid
    ngrid <- ngrid0 <- 10
    sl <- ceiling(ngrid*sqrt(maxgs.tot))
    ngrid <- ceiling(sqrt(sum(sl^2)))
    grid <- matrix(0,nrow=ngrid,ncol=ngrid) # generate initial grid
    gridsubs <- 1:(ngrid^2) # identifiers
    pointdata$xpos <- pointdata$ypos <- NA
    potpos <- c()
    blacklist <- c()
    maxcount <- 200

    ## extend grid
    for (g in 1:length(ord)) {
      subscripts <- pointdata$group==ord[g]
      count <- 0
      found <- FALSE
      while ((!found)&(count<=maxcount)) {
        if ((g>1)&(count==0)) {
          if (sl[ord[g]]<sl[ord[g-1]]) {
            potpos <- c()
            blacklist <- c()
          }
        }
        count2 <- 0
        while (length(potpos)==0) { # no potential position available
          if (count2>0) { # enlarge grid
            ngrid <- ngrid+1
            grid <- rbind(grid,rep(0,ngrid-1))
            grid <- cbind(grid,rep(0,ngrid))
            gridsubs <- seq(1,ngrid^2,1)
            ## correct blacklist
            blacklist <-  blacklist+floor(blacklist/ngrid)
          }
          ## determine potential grid positions
          potpos <- (grid[gridsubs]==0)& # occupied positions
          ((gridsubs%%ngrid)<(ngrid-sl[ord[g]]+2))& # to close bottom
          (gridsubs<(ngrid*(ngrid-sl[ord[g]]+1))) # to close left
          if (sl[ord[g]]>1) { # filter for groups with ns > 1 field
            potpos <- potpos&(gridsubs%%ngrid!=0)
          }
          potpos <- which(potpos)
          potpos <- potpos[!potpos%in%blacklist]
          count2 <- count2+1
        }

        xpotpos <- floor(potpos/ngrid)+1
        ypotpos <- potpos%%ngrid
        ypotpos[ypotpos==0] <- ngrid
        posind <- sample(1:length(potpos),1)
        pos <- potpos[posind]
        xpos <- xpotpos[posind]
        ypos <- ypotpos[posind]
        xsubs <- seq(xpos,xpos+sl[ord[g]]-1,1)
        ysubs <- seq(ypos,ypos+sl[ord[g]]-1,1)
        if (sum(c(grid[ysubs,xsubs])!=0)==0) {
          pointdata$xpos[subscripts] <- xpos
          pointdata$ypos[subscripts] <- ypos
          grid[ysubs,xsubs] <- ord[g]
          potpos <- potpos[!potpos%in%which(grid==ord[g])]
          found <- TRUE
        } else {
          count <- count+1
          blacklist <- c(blacklist,pos)
          potpos <- potpos[-posind]
        }
        if (count==maxcount) {
          warning(paste(" [!] found no grid position for group ",g,": ",as.character(traj.seqe)[group.id[group.name==2]],", omit!",sep=""))
        }
      }
    }
    pointdata$xpos <- (pointdata$xpos-1)/ngrid
    pointdata$ypos <- (pointdata$ypos-1)/ngrid

    ## determine widths of segments to draw

    pointdata[,paste("width",1:ncov,sep=".")] <-
      as.data.frame(sqrt(scale(x=pointdata[,paste("n",1:ncov,sep=".")],
                               center=FALSE,
                               scale=nid.cov.tot))*sqrt(grid.scale)*
                    cex*ngrid0/ngrid)

    ## determine central plot coordinates

    pointdata$xpos <- pointdata$xpos*sqrt(grid.scale)+
      1/2*(1-sqrt(grid.scale))+
        c(sqrt(maxgs.tot)*sqrt(grid.scale)*
          cex*ngrid0/ngrid)[pointdata$group]/2
    pointdata$ypos <- pointdata$ypos*sqrt(grid.scale)+
      1/2*(1-sqrt(grid.scale))+pointdata$y-0.5+
        c(sqrt(maxgs.tot)*sqrt(grid.scale)*
              cex*ngrid0/ngrid)[pointdata$group]/2

    ## ===================
    ## define group colors
    ## ===================

    if (is.null(cpal)) { # colors are not specified by the user
      cpal <-
        brewer.pal(8,"Dark2")[rep(seq(1,8),ceiling(ngroup/8))][1:ngroup]
    } else {
      cpal <- rep(cpal,length.out=ngroup)
    }
    cpal <- cpal[order(ord)]
    if (is.character(cpal)) {
      cpal <- col2rgb(cpal)
      cpal <- rgb(cpal[1,],cpal[2,],cpal[3,],
                  alpha*255,maxColorValue=255)
    }

    pointdata[,paste("col",1:ncov,sep=".")] <-
      matrix(rep(cpal[pointdata$group],ncov),ncol=ncov)
    if (!is.null(linedata)) {
      linedata[,paste("col",1:ncov,sep=".")] <-
        matrix(rep(cpal[linedata$group],ncov),ncol=ncov)
    }

    ## =================
    ## define line style
    ## =================

    if (!is.null(linedata)) {
      linedata$lty <- rep(lty,length.out=ngroup)[linedata$group]
    }

    ## ================
    ## execute the plot
    ## ================

    if (print) cat(" [>] plot the data\n")

    if (is.null(split)) {

      ## option no window split
      ## ======================

      splity <- 1 # dummy
      nsplits <- 1
    } else {

      ## window split
      ## ============

      ## option "by a certain order"
      if (is.numeric(split)) {
        pointdata$x <- pointdata$x-split
        subscripts.split <- pointdata$x==0
      }
      ## option "by first observation"
      if (split=="first") {
        subscripts.split <- pointdata$x==1
      }
      ## option "by last observation"
      if (split=="last") {
        for (g in unique(pointdata$group)) {
          pointdata$x[pointdata$group==g] <-
            pointdata$x[pointdata$group==g]-
              max(pointdata$x[pointdata$group==g])
          if (!is.null(linedata)) {
            if (sum(linedata$group==g)>0) {
              linedata[linedata$group==g,c("x0","x1")] <-
                linedata[linedata$group==g,c("x0","x1")]-
                  max(linedata$x1[linedata$group==g])
            }
          }
          xtickat <- xticklab <- min(pointdata$x):max(pointdata$x)
        }
        subscripts.split <- pointdata$x==0
      }
      splity <- sort(unique(pointdata$y[subscripts.split]))
    }

    if (!is.null(split)|!is.null(cov)) {
      nsplits <- ncov*length(splity)
      ## determine the plot layout
      if (is.null(layout)) {
        optpr <- optimpanelraster(nx=nsplits,
                                  ny=nsplits,
                                  npanels=nsplits,c=1)
        nxl <- optpr[1]
        nyl <- optpr[2]
      } else {
        nxl <- layout[2]
        nyl <- layout[1]
      }
      if (!((nxl==1)&(nyl==1))) {
        layout(matrix(1:(nxl*nyl),ncol=nxl,nrow=nyl,byrow=TRUE))
      }
    }

    ## set equal line widths in case of lweights are not demanded
    if (!lweights) {
      pointdata[,paste("width",1:ncov,sep=".")] <-
        min(pointdata[,paste("width",1:ncov,sep=".")])
    }

    ## set definitiv plot coordinates
    pointdata$xpos <- pointdata$xpos+pointdata$x-0.5
    if (!is.null(linedata)) {
      linedata <- merge(x=linedata,
                        y=pointdata[,c("group","x","y","xpos","ypos")],
                        by.x=c("group","x0","y0"),
                        by.y=c("group","x","y"),
                        all.x=TRUE,all.y=FALSE,sort=FALSE)
      linedata <- merge(x=linedata,
                        y=pointdata[,c("group","x","y","xpos","ypos")],
                        by.x=c("group","x1","y1"),
                        by.y=c("group","x","y"),
                        all.x=TRUE,all.y=FALSE,sort=FALSE)
      colnames(linedata)[(ncol(linedata)-3):ncol(linedata)] <-
        c("x0pos","y0pos","x1pos","y1pos")
    }

    ## set not observed trajectories to line width 0
    pointdata[,paste("width",1:ncov,sep=".")] <-
      apply(as.matrix(pointdata[,paste("width",1:ncov,sep=".")],ncol=ncov),
            2,zero2NA)
    if (!is.null(linedata)) {
      linedata[,paste("lwd",1:ncov,sep=".")] <-
        apply(as.matrix(linedata[,paste("lwd",1:ncov,sep=".")],ncol=ncov),
              2,zero2NA)
    }

    ## ======================
    ## show/hide trajectories
    ## ======================

    if (sum(!show%in%c(0,1))>0|show.type==3) {
      if (show.type==1) { # hide whole non-embeddable trajectories

        weights.group <- tapply(X=weights.use[traj.id]*w[traj.id],
                                INDEX=list(traj.group,cov.int[traj.id]),
                                FUN=sum)
        weights.group[is.na(weights.group)] <- 0
        weights.group <- scale(x=weights.group,center=FALSE,
                               scale=nid.cov.tot)
        hidem <- !((weights.group>=min(show))&(weights.group<=max(show)))
        pointdata[,paste("hide",1:ncov,sep=".")] <- hidem[pointdata$group,]
        if (!is.null(linedata)) {
          linedata[,paste("hide",1:ncov,sep=".")] <- hidem[linedata$group,]
        }
      } else if (show.type==2) { # hide parts of trajectories
        ## find indices to hide
        nrel.pd <-
          scale(as.matrix(pointdata[,paste("n",1:ncov,sep=".")],
                          ncol=ncov),center=FALSE,scale=nid.cov.tot)
        if (!is.null(linedata)) {
          nrel.ld <-
            scale(as.matrix(linedata[,paste("n",1:ncov,sep=".")],
                            ncol=ncov),center=FALSE,scale=nid.cov.tot)
          pointdata[,paste("hide",1:ncov,sep=".")] <-
            !((nrel.pd>=min(show))&(nrel.pd<=max(show)))
          linedata[,paste("hide",1:ncov,sep=".")] <-
            !((nrel.ld>=min(show))&(nrel.ld<=max(show)))
        }
      } else if (show.type==3) {
        if (ytype=="ordinal") {

          data.wide <- pointdata[,c("group","x","y")]
          data.wide$y <-
            factor(data.wide$y,levels=1:ny,labels=levels(y),ordered=TRUE)
          data.wide <-
            reshape(data.wide,timevar="x",idvar="group",direction="wide")
          data.wide <- data.wide[order(data.wide$group),]
          var.types <- list(ordratio=1:nx)
          dist <- as.matrix(daisy(data.wide[,-1],
                                  metric="gower",type=var.types))
          na.ins <- matrix(0,ncol=ncol(dist),nrow=nrow(dist))
          na.add <- apply(data.wide[,-1],1,function(x){sum(is.na(x))})
          na.ins <- na.ins+na.add
          na.ins[na.add==0,na.add!=0] <-
            na.ins[na.add==0,na.add!=0]+na.add[na.add!=0]
          dist <- dist+na.ins
          for (i in 1:ncov) {
            nid.group <-
              tapply(pointdata[,paste("n",i,sep=".")],pointdata$group,mean)
            m.nid.group <- nid.group%*%t(nid.group)
            diag(m.nid.group) <- diag(m.nid.group)-1
            m.nid.group[nid.group==0,nid.group==0] <- 0
            medoid <- rowSums(dist*m.nid.group)
            medoid[nid.group==0] <- Inf
            medoid <- which.min(medoid)
            data.wide[,paste("dmedoid",i,sep=".")] <- dist[medoid,]
            print(data.wide[medoid,])
          }
          data.wide[,paste("dmedoid",1:ncov,sep=".")] <-
            data.wide[,paste("dmedoid",1:ncov,sep=".")]/
              max(data.wide[,paste("dmedoid",1:ncov,sep=".")])
          hide.col <- col2rgb(hide.col)
          hide.col <- rgb2hsv(hide.col[1,],hide.col[2,],hide.col[3,])
          hide.col[3] <- 0
          hide.col <- hsv(hide.col[1,],hide.col[2,],hide.col[3,])
          pointdata[,paste("col",1:ncov,sep=".")] <- hide.col
          if (!is.null(linedata)) {
            linedata[,paste("col",1:ncov,sep=".")] <- hide.col
          }
          for (i in 1:ncov) {
            d.tmp <- data.wide[pointdata$group,paste("dmedoid",i,sep=".")]
            col.tmp <- col2rgb(pointdata[,paste("col",i,sep=".")])
            col.tmp <- rgb2hsv(col.tmp[1,],col.tmp[2,],col.tmp[3,])
            col.tmp[3,] <- d.tmp^hide.s
            col.tmp <- hsv(col.tmp[1,],col.tmp[2,],col.tmp[3,])
            pointdata[,paste("col",i,sep=".")] <- col.tmp
            if (!is.null(linedata)) {
              d.tmp <- data.wide[linedata$group,paste("dmedoid",i,sep=".")]
              col.tmp <- col2rgb(linedata[,paste("col",i,sep=".")])
              col.tmp <- rgb2hsv(col.tmp[1,],col.tmp[2,],col.tmp[3,])
              col.tmp[3,] <- d.tmp^hide.s
              col.tmp <- hsv(col.tmp[1,],col.tmp[2,],col.tmp[3,])
              linedata[,paste("col",i,sep=".")] <- col.tmp
            }
            maxgs <-
              as.matrix(1/data.wide[,paste("dmedoid",1:ncov,sep=".")])

          }
        } else {
          warning("show.type 3 is only available for ordinal data")
        }
      } else {
        pointdata[,paste("hide",1:ncov,sep=".")] <- FALSE
        if (!is.null(linedata)) {
          linedata[,paste("hide",1:ncov,sep=".")] <- FALSE
        }
      }
    }

    ## ======================
    ## define background grid
    ## ======================

    if (is.null(xlim)) {
      xlim <- c(min(pointdata$x)-sqrt(grid.scale)/2-
                sqrt(max((nid.cov.tot-nid.cov.use)/nid.cov.tot))*
                sqrt(grid.scale)*ngrid0/ngrid*cex,
                max(pointdata$x)+sqrt(grid.scale)/2)
    }

    if (is.null(ylim)) {
      ylim <- c(min(pointdata$y)-sqrt(grid.scale)/2-
               sqrt(max((nid.cov.tot-nid.cov.use)/nid.cov.tot))*
                sqrt(grid.scale)*ngrid0/ngrid*cex,
                max(pointdata$y)+sqrt(grid.scale)/2)
    }

    seggrid <- expand.grid(xgrid=
                           seq(min(c(ceiling(xlim[1]),min(pointdata$x))),
                               max(c(floor(xlim[2]),max(pointdata$x))),),
                           ygrid=1:ny)
    seggrid$col <- grid.fill

    for (sc in 1:ncov) {

      ## determine line orders
      ## =====================

      if ((show.type%in%c(1,2))&(sum(!show%in%c(0,1))>0)) {
        pg.hide <- which(tapply(!pointdata[,paste("hide",sc,sep=".")],
                                pointdata$group,sum)==0)
        plotgroups.cov <- setdiff(1:ngroup,pg.hide)
        plotgroups.cov <-
          plotgroups.cov[order(maxgs[plotgroups.cov,sc],
                               decreasing=lorder=="background")]
      } else {
        plotgroups.cov <-
          (1:ngroup)[order(maxgs[,sc],decreasing=lorder=="background")]
      }

      for (sy in splity) {

        ## extract groups with equal start events
        ## ======================================

        if (!is.null(split)) {
          plotgroups <-
            intersect(plotgroups.cov,
                      pointdata$group[subscripts.split&(pointdata$y==sy)])
        } else {
          plotgroups <- plotgroups.cov
        }

        ## call the plot
        ## =============

        plot(x=1,y=1,
             xlab=xlab,ylab=ylab,
             xlim=xlim,ylim=ylim,
             type="n",axes=FALSE,...)
        axis(1,xtickat,xticklab)
        axis(2,1:ny,levels(y),las=2)

        ## add titles
        ## ==========

        ## title for covariate
        if (mtext&(!is.null(cov)|!is.null(split)|
                   ((show.type%in%c(1,2))&
                    (sum(!show%in%c(0,1))>0)))) {
          mtextt <- NULL
          if (!is.null(cov)) {
            sub.sc <- paste("group = ",levels(cov)[sc],sep="")
            mtextt <- sub.sc
          }

          ## title text for category split
          if (!is.null(split)) {
            if (split%in%c("first","last")) {
              sub.sy <- paste(split," y = ",levels(y)[sy],sep="")
            }
            if (is.numeric(split)) {
              sub.sy <- paste(split,"=",levels(y)[sy])
            }
            mtextt <-
              paste(ifelse(is.null(mtextt),"",paste(mtextt,", ",sep="")),
                    sub.sy,sep="")
          }

          ## title for split proportion
          if ((show.type%in%c(1,2))&sum(!show%in%c(0,1))>0) {
            sub.show <- paste("rendered: ",
                              round(100*sum(maxgs[plotgroups,sc])+
                                    100*((nid.cov.tot-nid.cov.use)/
                                     nid.cov.tot)[sc],1),"%",sep="")

            mtextt <- paste(ifelse(is.null(mtextt),"",
                                  paste(mtextt,", ",sep="")),sub.show,sep="")
          }

          if (sy==1) {
            mtextt <- paste(mtextt,", n = ",nid.cov.tot[sc],sep="")
          }

          ## add title to plot
          if (!is.null(mtextt)) {
            if (((ncov*splity)>1)&(length(main)==1)) {
              title(main=mtextt)
            } else {
              mtext(text=mtextt,side=3)
            }
          }
        }

        if ((ncov*splity==1)) {
          title(main=main[1])
        } else {
          if (length(main)==1) {
            title(main=main[1],outer=TRUE)
          } else {
            print(main)
            if (length(main)>=(sc-1)*splity+sy) {
              title(main=main[(sc-1)*splity+sy],outer=FALSE)
            }
          }
        }

        ## draw a grid
        ## ===========

        if (grid.shape=="default") {
          rect(xleft=seggrid$xgrid-sqrt(grid.scale)/2,
               ybottom=seggrid$ygrid-sqrt(grid.scale)/2,
               xright=seggrid$xgrid+sqrt(grid.scale)/2,
               ytop=seggrid$ygrid+sqrt(grid.scale)/2,
               border=grid.fill,col=seggrid$col)
          if (grid.lwd!=0) {
            abline(h=1:ny,v=xtickat,col=grid.col,lwd=grid.lwd)
          }
        } else {
          if (grid.shape=="proportional") {
            rect(xleft=seggrid$xgrid-sqrt(grid.scale)/2,
                 ybottom=seggrid$ygrid-sqrt(grid.scale)/2,
                 xright=seggrid$xgrid+sqrt(grid.scale)/2,
                 ytop=seggrid$ygrid+sqrt(grid.scale)/2,
                 border=grid.fill,col=NULL)
            rect(xleft=seggrid$xgrid-sqrt(grid.scale)/2*ngrid0/ngrid*cex,
                 ybottom=seggrid$ygrid-sqrt(grid.scale)/2*ngrid0/ngrid*cex,
                 xright=seggrid$xgrid+sqrt(grid.scale)/2*ngrid0/ngrid*cex,
                 ytop=seggrid$ygrid+sqrt(grid.scale)/2*ngrid0/ngrid*cex,
                 border=grid.fill,col=grid.fill)
          }
        }

        ## plot the subjects without observations
        ## ======================================

        if (sy==1) {
          rect(xleft=min(pointdata$x)-sqrt(grid.scale)/2-
               sqrt((nid.cov.tot[sc]-nid.cov.use[sc])/nid.cov.tot[sc])*
               sqrt(grid.scale)*ngrid0/ngrid*cex,
               ybottom=min(pointdata$y)-sqrt(grid.scale)/2-
               sqrt((nid.cov.tot[sc]-nid.cov.use[sc])/nid.cov.tot[sc])*
               sqrt(grid.scale)*ngrid0/ngrid*cex,
               xright=min(pointdata$x)-sqrt(grid.scale)/2,
               ytop=min(pointdata$y)-sqrt(grid.scale)/2,col="black",
               border=0,lwd=0)
        }

        ## plot the trajectories
        ## =====================

        ## plot the hidden trajectories
        if ((show.type%in%c(1,2))&(sum(!show%in%c(0,1))>0)&
            (!is.null(hide.col))) {
          ss.p.hide <- c(!is.na(pointdata[,paste("width",sc,sep=".")]))&
          c(pointdata[,paste("hide",sc,sep=".")])
          rect(xleft=pointdata$xpos[ss.p.hide]-
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,,
               ybottom=pointdata$ypos[ss.p.hide]-
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,,
               xright=pointdata$xpos[ss.p.hide]+
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,
               ytop=pointdata$ypos[ss.p.hide]+
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,
               col=hide.col,border=border,lwd=border.lwd)
          ss.l.hide <- c(!is.na(linedata[,paste("lwd",sc,sep=".")]))&
          c(linedata[,paste("hide",sc,sep=".")])
          segments(x0=linedata$x0pos[ss.l.hide],
                   y0=linedata$y0pos[ss.l.hide],
                   x1=linedata$x1pos[ss.l.hide],
                   y1=linedata$y1pos[ss.l.hide],
                   lwd=linedata[ss.l.hide,paste("lwd",sc,sep=".")]*
                   (lwd.max-lwd.min)+lwd.min,col=hide.col,
                   lty=linedata$lty[ss.l.hide],lend=0)
        }

        ## plot the jittered trajectories
        if (type=="jitter") {
          ss.p <- (pointdata$group%in%plotgroups)&
          c(!is.na(pointdata[,paste("width",sc,sep=".")]))
          ss.d <- ss.p&(pointdata$Freq>1)
          if (!is.null(linedata)) {
            ss.l <- (linedata$group%in%plotgroups)&
            c(!is.na(linedata[,paste("lwd",sc,sep=".")]))
          }
          if (sum(ss.p)>0) {
            rect(xleft=pointdata$xpos[ss.p]-
                 c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                 ybottom=pointdata$ypos[ss.p]-
                 c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                 xright=pointdata$xpos[ss.p]+
                 c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                 ytop=pointdata$ypos[ss.p]+
                 c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                 col=pointdata[ss.p,paste("col",sc,sep=".")],
                 border=border,lwd=border.lwd)
          }
          if (!is.null(linedata)) {
            if (sum(ss.l)>0) {
              segments(x0=linedata$x0pos[ss.l],
                       y0=linedata$y0pos[ss.l],
                       x1=linedata$x1pos[ss.l],
                       y1=linedata$y1pos[ss.l],
                       lwd=linedata[ss.l,paste("lwd",sc,sep=".")]*
                       (lwd.max-lwd.min)+lwd.min,
                       col=linedata[ss.l,paste("col",sc,sep=".")],
                       lty=linedata$lty[ss.l],lend=0)
            }
          }
          if (sum(ss.d)>0) {
            points(x=pointdata$xpos[ss.d],
                   y=pointdata$ypos[ss.d],
                   pch=sf.pch,cex=sf.cex)
            i.multi <- which(ss.d)
            ppin <- par("pin")
            pusr <- par("usr")
            xr <- pointdata[ss.d,paste("width",sc,sep=".")]/2*
              sf.cex.leaves*abs(pusr[2L] - pusr[1L])/ppin[1L]
            yr <- pointdata[ss.d,paste("width",sc,sep=".")]/2*
              sf.cex.leaves*abs(pusr[4L] - pusr[3L])/ppin[2L]
            i.rep <- rep.int(i.multi,pointdata$Freq[ss.d])
            z <- numeric()
            for (i in i.multi) {
              z <- c(z, 1:pointdata$Freq[i] +
                     if (rotate) stats::runif(1) else 0)
            }
            deg <- (2 * pi * z)/pointdata$Freq[i.rep]

            segments(x0=pointdata$xpos[i.rep],
                     y0=pointdata$ypos[i.rep],
                     x1=pointdata$xpos[i.rep]+
                     xr*sin(deg),
                     y1=pointdata$ypos[i.rep]+
                     yr*cos(deg))
          }
        } else {
          ## plot the trajectories
          for (pg in plotgroups) {
            ss.p <- (pointdata$group%in%pg)&
            c(!is.na(pointdata[,paste("width",sc,sep=".")]))
            ss.d <- ss.p&(pointdata$Freq>1)
            if (!is.null(linedata)) {
              ss.l <- (linedata$group%in%pg)&
              c(!is.na(linedata[,paste("lwd",sc,sep=".")]))
            }

            ## draw the points
            if (sum(ss.p)>0) {
              rect(xleft=pointdata$xpos[ss.p]-
                   c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                   ybottom=pointdata$ypos[ss.p]-
                   c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                   xright=pointdata$xpos[ss.p]+
                   c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                   ytop=pointdata$ypos[ss.p]+
                   c(pointdata[ss.p,paste("width",sc,sep=".")])/2,
                   col=pointdata[ss.p,paste("col",sc,sep=".")],
                   border=border,lwd=border.lwd)
            }

            ## draw the lines
            if (!is.null(linedata)) {
              if (sum(ss.l)>0) {
                segments(x0=linedata$x0pos[ss.l],
                         y0=linedata$y0pos[ss.l],
                         x1=linedata$x1pos[ss.l],
                         y1=linedata$y1pos[ss.l],
                         lwd=linedata[ss.l,paste("lwd",sc,sep=".")]*
                         (lwd.max-lwd.min)+lwd.min,
                         col=linedata[ss.l,paste("col",sc,sep=".")],
                         lty=linedata$lty[ss.l],lend=0)
              }
            }

            ## inscribe doubles
            if (sum(ss.d)>0) {
              points(x=pointdata$xpos[ss.d],
                     y=pointdata$ypos[ss.d],
                     pch=sf.pch,cex=sf.cex)
              i.multi <- which(ss.d)
              ppin <- par("pin")
              pusr <- par("usr")
              xr <- pointdata[ss.d,paste("width",sc,sep=".")]/2*
                sf.cex.leaves*abs(pusr[2L] - pusr[1L])/ppin[1L]
              yr <- pointdata[ss.d,paste("width",sc,sep=".")]/2*
                sf.cex.leaves*abs(pusr[4L] - pusr[3L])/ppin[2L]
              i.rep <- rep.int(i.multi,pointdata$Freq[ss.d])
              z <- numeric()
              for (i in i.multi) {
                z <- c(z, 1:pointdata$Freq[i] +
                     if (rotate) stats::runif(1) else 0)
            }
              deg <- (2 * pi * z)/pointdata$Freq[i.rep]

              segments(x0=pointdata$xpos[i.rep],
                       y0=pointdata$ypos[i.rep],
                       x1=pointdata$xpos[i.rep]+
                       xr*sin(deg),
                       y1=pointdata$ypos[i.rep]+
                       yr*cos(deg))
            }
          }
        }
        box()
      }
    }

    if (return.data) {
      tmp <- pointdata[,c("Freq","x","y","group")]
      tmp <- data.frame(matrix(unlist(apply(tmp,1,reconstdoubles)),
                               ncol=3,byrow=TRUE))
      colnames(tmp) <- c("x","y","group")
      tmp <- tmp[order(tmp$group,tmp$x,tmp$y),]
      tmp$y <- factor(tmp$y,levels=1:ny,labels=levels(y))
      traj.string <- as.character(seqecreate(id=tmp$group,
                                             timestamp=tmp$x,event=tmp$y))
      size <- apply(X=as.matrix(pointdata[,c(paste("n",1:ncov,sep="."))],
                      ncol=ncov),
                    MARGIN=2,FUN=tapply,list(pointdata$group),max)
      if (!is.null(cov)) {
        dn2 <- paste("n",levels(cov),sep=".")
      } else {
        dn2 <- paste("n",1:ncov,sep=".")
      }
      ret <-
        matrix(rbind(nid.cov.tot-colSums(size),
                     matrix(size,ncol=ncov)),
               dimnames=list(c("()",traj.string),dn2),ncol=ncov)
      return(ret)
    }
  }

## ==================
## required functions
## ==================

## function for constructing time order
## ====================================

retorders <- function(x) {
  lev <- sort(unique(x),decreasing=FALSE)
  as.integer(factor(x,levels=lev))
}

## function for order alignement
## =============================

remfolltrues <- function(x) {
  ## erase all TRUE values following the forst TRUE value in
  ## a vector
  if (sum(x)>1) {
    x[which(x)[-1]] <- FALSE
  }
  return(x)
}

## erase 0 values for line widths and symbol sizes
## ===============================================

zero2NA <- function(x) {
  x[x==0] <- NA
  return(x)
}


## find optimal plot layout
## ========================

## optimization function

fmin <- function(nx,ny,n,c=0.5) {
  control <- nx*ny-n # number of empty windows
  opt <- control/n+c*abs(1-min(nx,ny)/max(nx,ny)) # optimization value:
                                        # function of number of empty windows
                                        # and ncol/nrow ratio
  return(c(control,opt))
}

## optimizer
optimpanelraster <- function(nx,ny,npanels,c=1) {
  minvalue1 <- fmin(nx-1,ny,npanels,c)
  minvalue2 <- fmin(nx,ny-1,npanels,c)
  while (minvalue1[1]>=0|minvalue2[1]>=0) {
    if ((minvalue1[1]>=0)&(minvalue2[1]>=0)) {
      if (minvalue1[2]<minvalue2[2]) {
        nx <- nx-1 } else { ny <- ny-1 }
    } else {
      if (minvalue1[2]>=0) {
        nx <- nx-1
      } else {
        if (minvalue2[2]>=0) {
          ny <- ny-1
        }
      }}
    minvalue1 <- fmin(nx-1,ny,npanels,c)
    minvalue2 <- fmin(nx,ny-1,npanels,c)
  }
  return(c(nx,ny))
}

## create trajectory strings
## =========================

createstring <- function(y,x=1:length(y),pa.left="(",pa.right=")",
                         pot="^",con="-") {
  subscripts <- !is.na(y);
  ret <- paste(pa.left,y,pa.right,pot,x,sep="");
  ret <- ret[subscripts];
  ret <- paste(ret,collapse=con)
  return(ret)
}

## reconstitute coordinate entries of simultaneous
## equal observations int the plot coordinates
## ===============================================

reconstdoubles <- function(x) {matrix(rep(x[2:length(x)],x[1]),nrow=x[1])}
