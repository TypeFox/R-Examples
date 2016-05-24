## ============================================
## convert between event sequence data formats
## ============================================

seqeformat <- function(data,from="TSE",to="seqe",
                       id=NULL,timestamp=NULL,event=NULL,
                       var=NULL,start=1,alphabet=NULL,states=NULL,
                       labels=NULL,weighted=TRUE,
                       weights=NULL,tevent="transition",
                       obs.intervals=NULL)
  {
    ## ==============================
    ## check obs.interval definition
    ## ==============================
    
    if (!is.null(obs.intervals)) {
      if (id!="id"&is.character(id)) {
        id <- which(colnames(obs.intervals)==id)
        colnames(obs.intervals)[id] <- "id"
      }
      if (sum(c("id","time.start","time.end")%in%
              colnames(obs.intervals))!=3) {
        warning(" [!] Invalid observation interval declaration. Ignore.")
        obs.intervals <- NULL
      }
    }

    ## ================
    ## treat TSE input
    ## ================ 

    ## prepare TSE input
    if (from=="TSE") {
      if (id!="id"&is.character(id)) {
        id <- which(colnames(data)==id)
        colnames(data)[id] <- "id"
      }
      data$id <- factor(data$id)
      if (timestamp!="timestamp"&is.character(timestamp)) {
        timestamp <- which(colnames(data)==timestamp)
        colnames(data)[timestamp] <- "timestamp"
      }
      if (event!="event"&is.character(event)) {
        event <- which(colnames(data)==event)
        colnames(data)[event] <- "event"
      }
      data$event <- factor(data$event)
    }

    ## convert TSE to TSE
    if ((from=="TSE")&(to=="TSE")) {
        if (!is.factor(data$id)) {
          data$id <- factor(data$id) }
        if (!is.factor(data$event)) {
          data$event <- factor(data$event) }
        data <- data[order(data$id,data$timestamp,data$event),]
        rownames(data) <- 1:nrow(data)
        ret <- data
      }

    ## convert TSE to seqe using seqecreate
    if ((from=="TSE")&(to=="seqe")) {
        ret <- seqecreate(data=data,id=id,timestamp=timestamp,
                          event=event,weighted=weighted)
      }

    ## convert TSE to both format (TSE and seqe)
    if ((from=="TSE")&(to=="both")) {
        if (!is.factor(data$id)) {
          data$id <- factor(data$id) }
        if (!is.factor(data$event)) {
          data$event <- factor(data$event) }
        data <- data[order(data$id,data$timestamp,data$event),]
        rownames(data) <- 1:nrow(data)
        ret <- list()
        ret$TSE <- data
        ret$seqe <- seqecreate(data=data,id=id,timestamp=timestamp,
                               event=event,weighted=weighted)       
      }

    ## ================
    ## treat STS input
    ## ================

    ## prepare STS input
    ## =================
    
    if (from=="STS") {
        if (is.null(id)) {seqid <- 1:nrow(data)} else {seqid <- data[,id]}
        data <- data[order(seqid),]
        seqid <- sort(seqid)
        ## data are STS raw data
        if (!inherits(data,"stslist")) {
          if (is.null(states)&!is.null(labels)) {
            states <- LETTERS[1:length(labels)]
            seq <- seqdef(data=data,var=var,informat="STS",
                          alphabet=alphabet,states=states,
                          id=seqid,labels=labels,weights=weights,
                          start=start)
          }
          if (is.null(states)) { # workaround due bug in seqdef
            seq <- seqdef(data=data,var=var,informat="STS",
                          alphabet=alphabet,
                          id=seqid,labels=labels,weights=weights,
                          start=start)
          }
          if (!exists("seq")) { # workaround due bug in seqdef
            seq <- seqdef(data=data,var=var,informat="STS",
                          alphabet=alphabet,states=states,
                          id=seqid,labels=labels,weights=weights,
                          start=start)
          }
        } else {
          ## data are formatted state sequencea
          seq <- data
        }
      }

    ## convert STS to seqe using seqcreate
    ## ===================================
    
    if ((from=="STS")&(to=="seqe")) {
      ret <- seqecreate(data=seq,weighted=weighted,
                        tevent=tevent)
    }

    ## convert STS to TSE using seqdef and seqformat
    ## =============================================
    
    if ((from=="STS")&(to=="TSE")) {
      ## create TSE data using seqformat
      m <- seqetm(seq=seq,method=tevent)
      ret <- seqformat(data=seq,from="STS",to="TSE",tevent=m)
      names(ret)[2] <- "timestamp"
      ret$timestamp <- ret$timestamp+start
      ret$id <- factor(rownames(seq)[ret$id],levels=rownames(seq))
      ## reorder event factor levels if period events were called
      if (tevent=="period") {
        which.end <- grep(pattern="end",x=levels(ret$event))
        which.start <- which(!1:nlevels(ret$event)%in%which.end)
        ind <- rep(c(1,1+length(which.start)),
                   length(which.start))+
                     rep(seq(0,length(which.start)-1),each=2)
        order <- c(which.start,which.end)[ind]           
        ret$event <- factor(ret$event,levels=levels(ret$event)[order])
      }
      ## add covariates from original data
      if (!inherits(data,"stslist")) {
        lockedvars <- which(colnames(data)%in%c("id","event","timestamp"))
        cov <- seq(1,ncol(data))[-var]
        if (length(cov)>1) {
          ret[,colnames(data)[cov]] <- data[ret$id,cov]
        }
      }
    }

    ## convert STS to both
    ## ===================
    
    if ((from=="STS")&(to=="both")) {
      ## create TSE data using seqformat
      ret <- list()
      m <- seqetm(seq=seq,method=tevent)
      ret$TSE <- seqformat(data=seq,from="STS",to="TSE",tevent=m)
      names(ret$TSE)[2] <- "timestamp"
      ret$TSE$timestamp <- ret$TSE$timestamp+start
      ret$TSE$id <- factor(rownames(seq)[ret$TSE$id],levels=rownames(seq))
      ret$seqe <- seqecreate(data=seq,weighted=weighted,tevent=tevent)
      ## reorder event factor levels if period events were called
      if (tevent=="period") {
        which.end <- grep(pattern="end",x=levels(ret$TSE$event))
        which.start <- which(!1:nlevels(ret$TSE$event)%in%which.end)
        ind <- rep(c(1,1+length(which.start)),
                   length(which.start))+
                     rep(seq(0,length(which.start)-1),each=2)
        order <- c(which.start,which.end)[ind]           
        ret$TSE$event <- factor(ret$TSE$event,
                                levels=levels(ret$TSE$event)[order])
      }
      ## add covariates from original data
      if (!inherits(data,"stslist")) {
        lockedvars <- which(colnames(data)%in%c("id","event","timestamp"))
        cov <- seq(1,ncol(data))[-c(var,lockedvars)]
        if (length(cov)>1) {
          ret$TSE[,colnames(data)[cov]] <- data[ret$TSE$id,cov]
        }
      }
      }
    
    ## ================================
    ## treat seqe input (experimental)
    ## ================================

    ## convert seqe to TSE or both
    if ((from=="seqe")&(to=="TSE")|
        (from=="seqe")&(to=="both"))
      {
        seqe.string <- as.character(data)
        split.seqe.string <- function(x){unlist(strsplit(x=x,split="-"))}
        seqe.decomp <- lapply(X=seqe.string,FUN=split.seqe.string)
        event <- vector("list",length(seqe.decomp))
        gaps <- vector("list",length(seqe.decomp))
        time <- vector("list",length(seqe.decomp))
        timerange <- vector("list",length(seqe.decomp))
        for (i in 1:length(seqe.decomp)) {
          if (substr(x=seqe.string[i],start=1,stop=1)=="(") {
            seqe.decomp[[i]] <- c("0",seqe.decomp[[i]])
          }
          event[[i]] <-
            seqe.decomp[[i]][seq(from=2,to=length(seqe.decomp[[i]]),by=2)]
          event[[i]] <- sub(pattern="(",replacement="",
                            x=event[[i]],fixed=TRUE)
          event[[i]] <- sub(pattern=")",replacement="",
                            x=event[[i]],fixed=TRUE)
          gaps[[i]] <-
            as.numeric(seqe.decomp[[i]][seq(from=1,
                                            to=length(seqe.decomp[[i]]),
                                            by=2)])
          timerange[[i]] <- c(gaps[[i]][1],sum(gaps[[i]]))
          if (length(gaps[[i]])>length(event[[i]]))
            {
              gaps[[i]] <- gaps[[i]][-length(gaps[[i]])]
            }
          time[[i]] <- cumsum(gaps[[i]])
          extract.simultanevents <- function(x){strsplit(x=x,split=",")}
          event[[i]] <- sapply(X=event[[i]],FUN=extract.simultanevents,
                               USE.NAMES=FALSE)
        }
        f.id <- function(x){length(unlist(x))}
        f.time <- function(x){unlist(lapply(X=x,FUN=length))}
        ret <- data.frame(id=factor(rep(1:length(event),
                            unlist(lapply(X=event,FUN=f.id)))),
                          timestamp=rep(unlist(time),
                            unlist(lapply(X=event,FUN=f.time))),
                          event=factor(unlist(event)))
      }
    
    ## convert seqe to both
    if ((from=="seqe")&(to=="both"))
      {
        ret <- list(TSE=ret,
                    seqe=data)
      }

    ## ==========================
    ## add observation time data
    ## ==========================
    
    if ((!is.null(obs.intervals))&(to=="both"))
      {
        obs.intervals <- obs.intervals[obs.intervals$id%in%
                                       levels(ret$TSE$id),]
        obs.intervals$id <- factor(obs.intervals$id,
                                   levels=levels(ret$TSE$id))
        ret$obs.intervals <- obs.intervals[order(obs.intervals$id),]
      }

    ## ===========================
    ## set class of output object
    ## ===========================
    
    if (to=="both")
      {
        class(ret) <- c("seqelist","list")
      }

    ## ==============
    ## return object
    ## ==============
    
    return(ret)
  }
