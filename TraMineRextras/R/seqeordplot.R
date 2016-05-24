## ===============================================
## Plotting
## ===============================================

seqeordplot <- function(seqe,group=NULL,weighted=TRUE,weights=NULL,
                        alphabet=NULL,
                        type="distinctive",embedding="most-frequent",
                        show=c(0,1),hide.col="grey75",
                        cpal=NULL,alpha=1,
                        lcourse="upwards",lorder="background",
                        lweights=TRUE,lwd.min=0.5,lwd.max=4,lty=1,
                        cex=1,border=grid.col,border.lwd=grid.lwd/2,
                        grid.col="white",grid.fill="grey90",
                        grid.scale=1/4,grid.shape="default",grid.lwd=0,
                        orderalign=ifelse(split=="last","last","first"),
                        split=NULL,layout=NULL,return.data=FALSE,
                        main="",sub=NULL,mtext=TRUE,
                        xlab="order position",ylab="",
                        xlim=NULL,ylim=NULL,...)


  {
    print <- FALSE

    ## ======
    ## checks
    ## ======

    if (!is.seqelist(seqe)) {
      stop(" [!] seqe should be a seqelist. See help on seqecreate.")
    }

    if (!(embedding%in%c("most-frequent","uniformly"))) {
      stop("[!] invalid specified type argument")
    }

    if (!lweights) {
      lwd$max <- lwd$min
    }

    if ((min(show)<0)|(max(show)>1)) {
      stop("[!] invalid specified show argmuent")
    }

    if (is.character(split)) {
      if (((orderalign=="last")&(split=="first"))|
          ((orderalign=="first")&(split=="last"))) {
        stop("[!] unallowed combination of orderalign and split arguments.")
      }
    }

    ## ==================
    ## preparing the data
    ## ==================

    if (print) cat(" [>] preparing the data\n")

    ## extract the event data
    ## ======================

    data <- seqe2OSE(seqe)
    id <- data$id
    x <- data$x
    y <- data$y
    rm(data)

    if (!is.null(group)) {
      if (length(group)!=nlevels(id)) {
        stop("[!] impossible to link id and group vector. See help")
      } else {
        names(group) <- levels(id)
      }
    }

    if (!is.null(alphabet)) {
      if (sum(levels(y) %in% alphabet) != nlevels(y)) {
        stop("[!] invalid alphabet.")
      } else {
        y <- factor(y,levels=alphabet)
      }
    }

    if (length(orderalign)>0) {
      if (!orderalign%in%c("first","last",levels(y))) {
        stop("[!] invalid specified orderalign argument")
      }
    }

    ## extract weights
    ## ===============

    if (weighted) {
      if (!is.null(weights)) {
        if (length(weights)!=nlevels(id)) {
          stop("[!] impossible to link id and weights vector. See help")
        } else {
          names(weights) <- levels(id)
        }
        if (length(unique(seqeweight(seqe)))!=1) {
          warning("[!] ignore weights stored in seqe object.")
        }
      } else {
        weights <- seqeweight(seqe)
      }
    } else {
      weights <- rep(1,nlevels(id))
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
        warning(paste(" [!] ignore NA entries of cases: ",
                      paste(subscripts.id,collapse=", "),
                      " in weights vector",sep=""))
        weights <- weights[!subscripts.ID]
        id <- factor(id[!subscripts.LONG],levels=names(weights))
        x <- x[!subscripts.LONG]
        y <- y[!subscripts.LONG]
        if (!is.null(group)) {
          group <- group[!subscripts.ID]
        }
      }
    }

    ## NAs in group variable
    if (!is.null(group)) {
      subscripts.ID <- is.na(group)
      subscripts.LONG <- id%in%names(group[subscripts.ID])
      if (sum(subscripts.LONG)>0) {
        subscripts.id <- setdiff(names(group),names(group[!subscripts.ID]))
        warning(paste(" [!] ignore NA entries of ids: ",
                      paste(subscripts.id,collapse=", "),
                      "in group variable",sep=""))
        group <- group[!subscripts.ID]
        weights <- weights[!subscripts.ID]
        id <- factor(id[!subscripts.LONG],levels=names(group))
        x <- x[!subscripts.LONG]
        y <- y[!subscripts.LONG]
      }
    }

    ## standardize weights

    if (!is.null(group)) {
      weights.group <- table(group)[as.integer(group)]
      weights.sum <-
        tapply(X=weights,INDEX=group,FUN=sum)[as.integer(group)]
      weights <- weights.use <- weights/weights.sum*weights.group
      names(weights) <- levels(id)
    } else {
      weights <- weights.use <- weights/sum(weights)*nlevels(id)
    }

    ## check the x variable
    ## ====================

    xtickat <- xticklab <- 1:max(x)

    ## check the id variable
    ## =====================

    nid.tot <- nid.use <- nlevels(id)
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

    ## check group variable
    ## ====================

    if (!is.null(group)) {
      if (!is.factor(group)) {
        group <- factor(group)
      }
      nid.group.tot <- nid.group.use <- table(group)
      if (sum(nid.group.tot==0)>0) {
        warning(paste(" [!] the plot will omit the following group categories:",paste(names(nid.group.tot)[nid.group.tot==0],collapse=", "),sep=" "))
        group <- factor(group,levels=names(nid.group.tot)[nid.group.tot>0])
        nid.group.tot <- nid.group.use <- table(group)
      }
      ngroup <- nlevels(group)
      group.int <- as.integer(group)
      group.old <- group
    } else {
      ngroup <- 1
      group.int <- rep(1,nid.use)
      nid.group.tot <- nid.tot
      nid.group.use <- nid.use
    }

    ## order alignment
    ## ===============

    if (is.character(orderalign)) {
      if (orderalign=="last") {
        align <- tapply(X=x,INDEX=id.int,FUN=max)
        x <- x-align[id.int]
        x <- x-min(x)+1
        xtickat <- 1:max(x)
        xticklab <- seq(1-max(x),0)
      } else if (orderalign%in%levels(y)) {
        subscripts <- y==orderalign
        subscripts <- unlist(tapply(X=subscripts,
                                    INDEX=list(id.int),
                                    FUN=remfolltrues))
        aglevs <- c("no align",paste(orderalign,"align"))
        align <- rep(1,nid.use)
        tab <- table(y,id)[orderalign,]
        align[which(tab>0)] <- x[subscripts]
        aligngroup <- factor(aglevs[(tab>0)+1])
        if (!is.null(group)) {
          group <- factor(paste(group,aligngroup,sep="/"))
        } else {
          group <- aligngroup
        }
        ngroup <- nlevels(group)
        group.int <- as.integer(group)
        nid.group.tot <- nid.group.use <- table(group)
        x <- x-align[id.int]
        x <- x-min(x)+1
        xtickat <- 1:max(x)
        xticklab <- xtickat-x[which(y==orderalign)[1]]
      }
    }

    nx <- max(x)

    ## construct x and y lists (one element per id)
    ## ============================================

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
    ## in case a groupariate vector was entered
    ## ==================================================

    tab <- table(traj.id,traj.group)
    if ((type=="non-embeddable")&(ngroup>1)) {
      uemb <- rowSums(tab)==1
      for (i in 1:ngroup) {
        ## for which non-embeddedable trajectories there are observed
        ## trajectories that were solely embedded into it

        uemb.groups <- traj.group[traj.id%in%which(uemb&(group.int==i))]
        for (j in which((group.int==i)&(!uemb)))
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

    ## extract coordinate pairs for plotting the lines using
    ## segments (necessary to allow different line widths)
    ## =====================================================

    ntraj <- max(pointdata$group)
    linedata <- NULL

    for (g in 1:ntraj) {
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

    pointdata[,paste("n",1:ngroup,sep=".")] <- 1
    if (!is.null(linedata)) {
      linedata[,paste("n",1:ngroup,sep=".")] <- 0
      linedata[,paste("lwd",1:ngroup,sep=".")] <- 1
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
                                factor(group.int[w==1],levels=1:ngroup)),
                              FUN=sum)
      weights.group[is.na(weights.group)] <- 0

      ## perform unique group assignements
      for (i in which(w!=1)) {
        group.iin <- traj.group[traj.id==i]
        group.iin <-
          group.iin[which.max(weights.group[group.iin,group.int[i]])]
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
      pointdata[i,paste("n",1:ngroup,sep=".")] <-
        unlist(tapply(X=weights.use[isubs]*w[isubs]*
                      ((minx[isubs]<=pointdata$x[i])&
                       (maxx[isubs]>=pointdata$x[i])),
                      INDEX=list(factor(group.int[isubs],levels=1:ngroup)),
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
        linedata[i,paste("n",1:ngroup,sep=".")] <-
          unlist(tapply(X=weights.use[isubs]*w[isubs]*
                        ((minx[isubs]<=linedata$x0[i])&
                         (maxx[isubs]>=linedata$x1[i])),
                        INDEX=list(factor(group.int[isubs],levels=1:ngroup)),
                        FUN=sum))
      }
      linedata[is.na(linedata)] <- 0
    }

    ## line widths
    if (!is.null(linedata)) {
      if (!is.null(group)) {
        maxprop <-
          apply(as.matrix(linedata[,paste("n",1:ngroup,sep=".")],
                          ncol=ngroup),2,max,na.rm=TRUE)/nid.group.tot

        facm <- diag(maxprop/max(maxprop))
        scalem <- apply(as.matrix(linedata[,paste("n",1:ngroup,sep=".")],
                                  ncol=ngroup),2,max,na.rm=TRUE)
        scalem[scalem==0] <- 1
        scalem <- scale(x=linedata[,paste("n",1:ngroup,sep=".")],
                        center=FALSE,scale=scalem)
        scalem <- scalem%*%facm
        linedata[,paste("lwd",1:ngroup,sep=".")] <- as.data.frame(scalem)
      } else {
        linedata[,paste("lwd",1:ngroup,sep=".")] <-
          linedata[,paste("n",1:ngroup,sep=".")]/
            max(linedata[,paste("n",1:ngroup,sep=".")])
      }
    }

    ## ==============================
    ## Find plot positions for groups
    ## ==============================

    if (print) cat(" [>] determine plot coordinates\n")

    maxgs <- apply(X=as.matrix(x=pointdata[,paste("n",1:ngroup,sep=".")],
                     ncol=ngroup),MARGIN=2,FUN=tapply,
                   list(pointdata$group),max)
    maxgs <- scale(x=matrix(maxgs,ncol=ngroup),
                   center=FALSE,scale=nid.group.tot)
    maxgs.tot <- apply(X=as.matrix(x=maxgs,ncol=ngroup),MARGIN=1,FUN=max)
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

    pointdata[,paste("width",1:ngroup,sep=".")] <-
      as.data.frame(sqrt(scale(x=pointdata[,paste("n",1:ngroup,sep=".")],
                               center=FALSE,
                               scale=nid.group.tot))*sqrt(grid.scale)*
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
        brewer.pal(8,"Dark2")[rep(seq(1,8),ceiling(ntraj/8))][1:ntraj]
    } else {
      cpal <- rep(cpal,length.out=ntraj)
    }
    cpal <- cpal[order(ord)]
    if (is.character(cpal)) {
      cpal <- col2rgb(cpal)
      cpal <- rgb(cpal[1,],cpal[2,],cpal[3,],
                  alpha*255,maxColorValue=255)
    }

    pointdata[,paste("col",1:ngroup,sep=".")] <-
      matrix(rep(cpal[pointdata$group],ngroup),ncol=ngroup)
    if (!is.null(linedata)) {
      linedata[,paste("col",1:ngroup,sep=".")] <-
        matrix(rep(cpal[linedata$group],ngroup),ncol=ngroup)
    }

    ## =================
    ## define line style
    ## =================

    if (!is.null(linedata)) {
      linedata$lty <- rep(lty,length.out=ntraj)[linedata$group]
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
        if (!is.null(linedata)) {
          linedata[,c("x0","x1")] <- linedata[,c("x0","x1")]-split
        }
        xtickat <- min(pointdata$x):max(pointdata$x)
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

    if (!is.null(split)|!is.null(group)) {
      nsplits <- ngroup*length(splity)
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
      pointdata[,paste("width",1:ngroup,sep=".")] <-
        min(pointdata[,paste("width",1:ngroup,sep=".")])
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
    pointdata[,paste("width",1:ngroup,sep=".")] <-
      apply(as.matrix(pointdata[,paste("width",1:ngroup,sep=".")],
                      ncol=ngroup),2,zero2NA)
    if (!is.null(linedata)) {
      linedata[,paste("lwd",1:ngroup,sep=".")] <-
        apply(as.matrix(linedata[,paste("lwd",1:ngroup,sep=".")],
                        ncol=ngroup),2,zero2NA)
    }

    ## ======================
    ## show/hide trajectories
    ## ======================

    if (sum(!show%in%c(0,1))>0) {
        weights.group <- tapply(X=weights.use[traj.id]*w[traj.id],
                                INDEX=list(traj.group,group.int[traj.id]),
                                FUN=sum)
        weights.group[is.na(weights.group)] <- 0
        weights.group <- scale(x=weights.group,center=FALSE,
                               scale=nid.group.tot)
        hidem <- !((weights.group>=min(show))&(weights.group<=max(show)))
        pointdata[,paste("hide",1:ngroup,sep=".")] <- hidem[pointdata$group,]
        if (!is.null(linedata)) {
          linedata[,paste("hide",1:ngroup,sep=".")] <- hidem[linedata$group,]
        }
      }

    ## ======================
    ## define background grid
    ## ======================

    if (is.null(xlim)) {
      xlim <- c(min(pointdata$x)-sqrt(grid.scale)/2-
                sqrt(max((nid.group.tot-nid.group.use)/nid.group.tot))*
                sqrt(grid.scale)*ngrid0/ngrid*cex,
                max(pointdata$x)+sqrt(grid.scale)/2)
    }

    if (is.null(ylim)) {
      ylim <- c(min(pointdata$y)-sqrt(grid.scale)/2-
               sqrt(max((nid.group.tot-nid.group.use)/nid.group.tot))*
                sqrt(grid.scale)*ngrid0/ngrid*cex,
                max(pointdata$y)+sqrt(grid.scale)/2)
    }

    seggrid <- expand.grid(xgrid=
                           seq(min(c(ceiling(xlim[1]),min(pointdata$x))),
                               max(c(floor(xlim[2]),max(pointdata$x))),),
                           ygrid=1:ny)
    seggrid$col <- grid.fill

    for (sc in 1:ngroup) {

      ## determine line orders
      ## =====================

      if (sum(!show%in%c(0,1))>0) {
        pg.hide <- which(tapply(!pointdata[,paste("hide",sc,sep=".")],
                                pointdata$group,sum)==0)
        plotgroups.group <- setdiff(1:ntraj,pg.hide)
        plotgroups.group <-
          plotgroups.group[order(maxgs[plotgroups.group,sc],
                                 decreasing=lorder=="background")]
      } else {
        pg.hide <- c()
        plotgroups.group <-
          (1:ntraj)[order(maxgs[,sc],decreasing=lorder=="background")]
      }

      for (sy in splity) {

        ## extract groups with equal start events
        ## ======================================

        if (!is.null(split)) {
          plotgroups <-
            intersect(plotgroups.group,
                      pointdata$group[subscripts.split&(pointdata$y==sy)])
          hidegroups <-
            intersect(pg.hide,
                      pointdata$group[subscripts.split&(pointdata$y==sy)])
        } else {
          plotgroups <- plotgroups.group
          hidegroups <- pg.hide
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

        ## title for groupariate
        if (mtext&(!is.null(group)|!is.null(split)|
                   (sum(!show%in%c(0,1))>0))) {
          mtextt <- NULL
          if (!is.null(group)) {
            sub.sc <- paste("group = ",levels(group)[sc],sep="")
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
          if (sum(!show%in%c(0,1))>0) {
            sub.show <- paste("rendered: ",
                              round(100*sum(maxgs[plotgroups,sc])+
                                    100*((nid.group.tot-nid.group.use)/
                                     nid.group.tot)[sc],1),"%",sep="")

            mtextt <- paste(ifelse(is.null(mtextt),"",
                                  paste(mtextt,", ",sep="")),sub.show,sep="")
          }

          if (sy==splity[1]) {
            mtextt <- paste(mtextt,", n = ",nid.group.tot[sc],sep="")
          }

          ## add title to plot
          if (!is.null(mtextt)) {
            if (((ngroup*nsplits)>1)&(length(main)==1)) {
              title(main=mtextt)
            } else {
              mtext(text=mtextt,side=3)
            }
          }
        }

        if ((ngroup*nsplits==1)) {
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

        if (sy==splity[1]) {
          rect(xleft=min(pointdata$x)-sqrt(grid.scale)/2-
               sqrt((nid.group.tot[sc]-nid.group.use[sc])/
                    nid.group.tot[sc])*sqrt(grid.scale)*ngrid0/ngrid*cex,
               ybottom=min(pointdata$y)-sqrt(grid.scale)/2-
               sqrt((nid.group.tot[sc]-nid.group.use[sc])/
                    nid.group.tot[sc])*sqrt(grid.scale)*ngrid0/ngrid*cex,
               xright=min(pointdata$x)-sqrt(grid.scale)/2,
               ytop=min(pointdata$y)-sqrt(grid.scale)/2,col="black",
               border=0,lwd=0)
        }

        ## plot the trajectories
        ## =====================

        ## plot the hidden trajectories
        if ((sum(!show%in%c(0,1))>0)&(!is.null(hide.col))) {
          ss.p.hide <- (pointdata$group%in%hidegroups)&
          c(!is.na(pointdata[,paste("width",sc,sep=".")]))
          rect(xleft=pointdata$xpos[ss.p.hide]-
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,,
               ybottom=pointdata$ypos[ss.p.hide]-
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,,
               xright=pointdata$xpos[ss.p.hide]+
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,
               ytop=pointdata$ypos[ss.p.hide]+
               c(pointdata[ss.p.hide,paste("width",sc,sep=".")])/2,
               col=hide.col,border=border,lwd=border.lwd)
          ss.l.hide <- (linedata$group%in%hidegroups)&
          c(!is.na(linedata[,paste("lwd",sc,sep=".")]))
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
        } else {
          ## plot the trajectories
          for (pg in plotgroups) {
            ss.p <- (pointdata$group%in%pg)&
            c(!is.na(pointdata[,paste("width",sc,sep=".")]))
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
          }
        }
        box()
      }
    }

    if (return.data) {
      tmp <- pointdata[,c("x","y","group")]
      tmp <- tmp[order(tmp$group,tmp$x,tmp$y),]
      tmp$y <- factor(tmp$y,levels=1:ny,labels=levels(y))
      traj.string <- as.character(seqecreate(id=tmp$group,
                                             timestamp=tmp$x,event=tmp$y))
      size <- apply(X=as.matrix(pointdata[,c(paste("n",1:ngroup,sep="."))],
                      ncol=ngroup),
                    MARGIN=2,FUN=tapply,list(pointdata$group),max)
      if (!is.null(group)) {
        dn2 <- paste("n",levels(group),sep=".")
      } else {
        dn2 <- paste("n",1:ngroup,sep=".")
      }
      ret <-
        matrix(rbind(nid.group.tot-colSums(size),matrix(size,ncol=ngroup)),
               dimnames=list(c("()",traj.string),dn2),ncol=ngroup)
      return(ret)
    }
  }

## ==================
## required functions
## ==================

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

## convert seqe object to order stamped event format
## =================================================

seqe2OSE <- function(seqe) {
  seqe.string <- as.character(seqe)
  seqe.string <-
    unlist(lapply(X=seqe.string, # erase initial time points
                  FUN=function(x){
                    substr(x=x,start=grep("(",x,fixed=TRUE),
                           stop=nchar(x))}))
  seqe.decomp <- lapply(X=seqe.string,
                        FUN=function(x){unlist(strsplit(x=x,split="-"))})
  seqe.decomp <-
    lapply(X=seqe.decomp,FUN=function(x){x[seq(1,length(x),2)]})
  seqe.decomp <- lapply(X=seqe.decomp,
                        FUN=function(x){sub(pattern="(",replacement="",
                          x=x,fixed=TRUE)})
  seqe.decomp <- lapply(X=seqe.decomp,
                        FUN=function(x){sub(pattern=")",replacement="",
                          x=x,fixed=TRUE)})
  seqe.decomp <- lapply(X=seqe.decomp,
                        FUN=function(x){
                          sapply(X=x,
                                 FUN=function(x){strsplit(x=x,split=",")},
                                 USE.NAMES=FALSE)})
  order <-
    lapply(X=seqe.decomp,FUN=function(x){unlist(lapply(X=x,FUN=length))})
  order <- lapply(X=order,FUN=function(x){rep(seq(1,length(x)),x)})
  id <- factor(rep(seq(1:length(order)),unlist(lapply(X=order,FUN=length))))
  order <- unlist(order)
  event <- unlist(seqe.decomp)
  return(data.frame(id=id,x=order,y=event))
}
