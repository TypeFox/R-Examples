##
## MDR
## 
## by SungGON Yi. <skon@kr.FreeBSD.org> and Eun-kyung Lee
##
## $Id$
##

##dataset <- read.csv(file1)

## dataset object must be data.frame (or list).
## response data must be coded by 1 (case) or 0 (control).
## snp data must be coded by 0, 1, or 2.
.PACKAGE <- "ORMDR" 

mdr.c <- function(dataset, colresp, cs, combi, cv.fold = 10, randomize = TRUE) {

	errRate.C <- function(comb, train, test, threshold) {
	  z <- .C("err_rate", as.integer(comb), nrow(comb), ncol(comb),
	          as.integer(unlist(train)), nrow(train), ncol(train),
	          as.integer(unlist(test)), nrow(test),
	          as.double(threshold), 
	          err.train = double(ncol(comb)), 
	          err.test = double(ncol(comb)), NAOK = FALSE,PACKAGE="ORMDR")
	  cbind(train = z$err.train, test = z$err.test)
	}  

    resp <- dataset[, colresp]
    case <- which(resp == cs)
    ctl <- which(resp != cs)
    if (randomize) {
        case <- sample(case)
        ctl <- sample(ctl)
    }
    resp <- as.integer(resp == cs)
    snp <- dataset[, -colresp]

    cv.case <- matrix(0, nrow = length(case), ncol = 2)
    cv.case[, 1] <- 1:length(case) %% cv.fold + 1
#  cv.case[, 2] <- sample(case)
    cv.case[, 2] <- case
    cv.ctl <- matrix(0, nrow = length(ctl), ncol = 2)
    cv.ctl[, 1] <- 1:length(ctl) %% cv.fold + 1
#  cv.ctl[, 2] <- sample(ctl)
    cv.ctl[, 2] <- ctl
    cv <- rbind(cv.case, cv.ctl)
    d <- cbind(resp, snp)

    ## result
    z <- list()

    ## combinations
    comb <- list()
    k <-combi
    comb <- combn(ncol(snp), k)
    comb <- comb + 1
    z <- list()
    z$min.comb <- matrix(0, nrow = k, ncol = cv.fold)
    z$train.erate <- numeric(cv.fold)
    z$test.erate <- numeric(cv.fold)


    for (i in 1:cv.fold) {
        train <- d[-cv[cv[, 1] == i, 2], ]
        test <- d[cv[cv[, 1] == i, 2], ]
        threshold <- length(which(train[, 1] == cs)) / 
        length(which(train[, 1] != cs)) 

        ## train and test error ate for all combination of SNPs
        errs <- errRate.C(comb, train, test, threshold)
        min.loc <- which.min(errs[, 1])
        z$min.comb[, i] <- comb[, min.loc] - 1
        z$train.erate[i] <- errs[min.loc, 1]
        z$test.erate[i] <- errs[min.loc, 2]
    }
    z$min.comb[z$min.comb >= colresp] <- z$min.comb[z$min.comb >= colresp] + 1

    ## maximum repeated selection
  
    z$data <- d
 
    name<-apply(t(z$min.comb),1,function(x) {
        t<-NULL;
        for(i in 1:length(x)) 
            t<-paste(t,as.character(x[i]),sep=";") ;return(t)})
    z$best.combi<-unlist(strsplit(names(sort(table(name),decreasing=TRUE))[1],";"))[-1]

    z
}

ormdr<-function(dataset,bestcombi,cs,colresp,CI.Asy=TRUE,CI.Boot=FALSE,B=5000) {  
    case.id<-which(dataset[,colresp]==cs)
    control.id<-c(1:nrow(dataset))[-case.id]
    snp <- dataset[, -colresp]
    bestcombi0 <- bestcombi
    bestcombi0[bestcombi >= colresp] <- bestcombi[bestcombi >= colresp] - 1
  
    best.data.case<-list()
    best.data.cont<-list()
    for(i in 1:length(bestcombi)) {
#        best.data.case[[i]]<-factor(dataset[case.id,bestcombi[i]],levels=c(0,1,2))
#        best.data.cont[[i]]<-factor(dataset[control.id,bestcombi[i]],levels=c(0,1,2))
        best.data.case[[i]]<-factor(snp[case.id,bestcombi0[i]],levels=c(0,1,2))
        best.data.cont[[i]]<-factor(snp[control.id,bestcombi0[i]],levels=c(0,1,2))
    }
   
    n.case<-length(case.id)
    n.cont<-length(control.id)
    t.case<-table(best.data.case)/n.case
    t.control<-table(best.data.cont)/n.cont
    Odds<-c(t.case/t.control)
   
#    LU.Asy<-LU.Boot<-matrix(NA,ncol=2,nrow=3**length(bestcombi))
    if(CI.Asy) {
        LU.Asy <-matrix(NA,ncol=2,nrow=3**length(bestcombi))
        L<-c(exp(log(Odds)-1.96*sqrt((1-t.case)/(t.case*n.case)+(1-t.control)/(t.control*n.cont))))
        U<-c(exp(log(Odds)+1.96*sqrt((1-t.case)/(t.case*n.case)+(1-t.control)/(t.control*n.cont))))
        LU.Asy<-cbind(L,U)
    }
    if(CI.Boot) {
        LU.Boot<-matrix(NA,ncol=2,nrow=3**length(bestcombi))
        Odds.keep<-NULL
        for(i in 1:B) {
            case.boot<-sample(1:n.case,n.case,replace=TRUE)
            cont.boot<-sample(1:n.cont,n.cont,replace=TRUE)
            boot.data.case<-list()
            boot.data.cont<-list()
            for(i in 1:length(bestcombi)) {
#                boot.data.case[[i]]<-factor(dataset[case.boot,bestcombi[i]],levels=c(0,1,2))
#                boot.data.cont[[i]]<-factor(dataset[control.boot,bestcombi[i]],levels=c(0,1,2))
                boot.data.case[[i]]<-factor(snp[case.boot,bestcombi0[i]],levels=c(0,1,2))
                boot.data.cont[[i]]<-factor(snp[cont.boot,bestcombi0[i]],levels=c(0,1,2))
            }
            P1.t<-c(table(boot.data.case)/n.case)
            P2.t<-c(table(boot.data.cont)/n.cont)
            Odds.t<-P1.t/P2.t
            Odds.keep<-rbind(Odds.keep,c(Odds.t))
        }
        LU.Boot<-NULL
        for(id in 1:ncol(Odds.keep)) {
            L.t<-sort(Odds.keep[,id])[round(B*0.025)]
            U.t<-sort(Odds.keep[,id])[round(B*0.975)]
            LU.Boot<-rbind(LU.Boot,c(L.t,U.t))
        }

    }
#    classID<-cbind(rep(0:2,3),rep(0:2,each=3))
#    if(length(bestcombi)>2) {
#        for(i in 3:length(bestcombi))
#            classID<-cbind(rbind(classID,classID,classID),rep(0:2,each=nrow(classID)))
#    }
    classID <- sapply(1:length(bestcombi), function(x) gl(3, 3^(x - 1), 3^length(bestcombi), labels = 0:2))
#    storage.mode(classID) <- "integer"
    cell.freq<-cbind(c(table(best.data.case)),c(table(best.data.cont)))
    Hi.Low<-ifelse(Odds>=1,"High","Low")
#     ORMDR.table<-cbind(classID,cell.freq,Hi.Low,round(Odds,3),rank(Odds),round(LU.Asy,3),round(LU.Boot,3))
    ORMDR.table<-cbind(classID,cell.freq,Hi.Low,round(Odds,3),rank(Odds))
    ORMDR.name <- c(colnames(snp)[bestcombi0], "case.freq","cont.freq","Hi.Low","Odds.ratio","Rank")
    if (CI.Asy) {
        ORMDR.table <- cbind(ORMDR.table, round(LU.Asy,3))
        ORMDR.name <- c(ORMDR.name, "Asy.L", "Asy.U")
    }
    if (CI.Boot) {
        ORMDR.table <- cbind(ORMDR.table, round(LU.Boot,3))
        ORMDR.name <- c(ORMDR.name, "Boot.L", "Boot.U")
    }
     
#    colnames(ORMDR.table)<-c(colnames(dataset)[bestcombi],                                        "case.freq","cont.freq","Hi.Low","Odds.ratio","Rank","Asy.L","Asy.U","Boot.L","Boot.U")
    colnames(ORMDR.table) <- ORMDR.name
    
    ORMDR.table
}

 .onLoad<-function(lib,pkg){
  library.dynam("ORMDR",pkg,lib)
}
