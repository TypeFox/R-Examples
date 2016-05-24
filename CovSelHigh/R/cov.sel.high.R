cov.sel.high<-function(T=NULL, Y=NULL, X=NULL,type=c("mmpc","mmhc","rf","lasso"),  betahat=TRUE, parallel=FALSE, Simulate=TRUE, N=NULL, Setting=1, rep=1, Models=c("Linear", "Nonlinear", "Binary"),...){
  Simulate<<-substitute(Simulate)
  N<<-substitute(N)
  Setting<<-substitute(Setting)
  Rep<<-substitute(rep)
  Models<<-match.arg(Models)
  type<<-match.arg(type)
  betahat<<-substitute(betahat)


  if(Simulate==FALSE){
                        if(is.null(T)==FALSE && is.null(Y)==FALSE && is.null(X)==FALSE){
                              Setting<-Models<-NULL
                              Rep<-1
                              N<-dim(X)[1]
                              }else{stop("Data or simulation settings must be provided correctly")}

                        if (sum(is.na(X)) > 0 | sum(is.na(T)) > 0 | sum(is.na(Y)) > 0) {
                            stop("missing data is currently not supported. Check T, Y, and X for missing values")
                            }

                        if (class(Y)!="numeric") {
                            stop("the outcome Y must be a numeric vector")
                            }

                        if (length(unique(T))!=2 || is.na(match(unique(T),c(0,1))) ){
                            stop("the treatment variable T must be binary")
                            }
                      uniqueclass<-unique(unlist(lapply(X,class)))
                      nrclass<-length(uniqueclass)
                      if(sum(uniqueclass %in% c("factor", "ordered", "numeric"))<nrclass){
                              stop("the covariates in X must be either numeric, factors or ordered factors.")
                            }

                        dat<-data.frame(X,Y,T)
                      }

    if(Simulate==TRUE){
                        if(is.null(N) || is.null(Setting) || is.null(Rep) || is.null(Models)){stop("Data or simulation settings must be provided correctly")}
                        if(parallel==TRUE){registerDoRNG(1231255125)}else{set.seed(1231255125)}

                      }

    reslist<-vector("list", Rep)

    if(parallel==TRUE){

      reslist<-foreach(i=1:Rep, .packages=c('MASS', 'bindata','Matching','bnlearn','CovSelHigh','glmnet','randomForest'),.export=c("Simulate","N","Setting","Rep","Models","type","betahat"))%dopar%{
       
          if(Simulate==TRUE){
            dat<-cov.sel.high.sim(N, Setting, Rep, Models)$dat
            X<-dat[,-c((dim(dat)[2])-1,dim(dat)[2])]
            Y<-dat[,((dim(dat)[2])-1)]
            T<-dat[,dim(dat)[2]]
          }
          varnames<-colnames(dat)
          
          covarcol<-1:(dim(X)[2])
          ycol<-(dim(X)[2])+1
          Tcol<-ycol+1
          if(class(dat[,Tcol])!="factor"){
            dat[,Tcol]<-factor(dat[,Tcol])
          }
          
          datbeta<-dat
          if(length(unique(dat[,ycol]))==2){
            datbeta[,ycol]<-factor(datbeta[,ycol])
          }
          
          
          
          
          
          
          numcol<-which(lapply(dat[covarcol],class)=="numeric")
          
          
          covars<-colnames(dat[,covarcol])
          
          
          if(type=="mmpc" || type=="mmhc"){
            if(length(unique(dat[,ycol]))>2){
              tried<-try(discretize(data.frame(dat[,ycol]),method="quantile"),silent=TRUE)
              if(class(tried)=="data.frame"){
                dat[,ycol]<-discretize(data.frame(dat[,ycol]),method="quantile")
              }else{
                print(tried)
                stop("the numeric outcome could not be discretized, recode it into factor")
                
              }
            }else{dat[,ycol]<-factor(dat[,ycol])
            }
            
            if(length(numcol)>0){
              if(length(numcol)==1){
                tried<-try(discretize(data.frame(dat[,numcol]),method="quantile"),silent=TRUE)
                if(class(tried)=="data.frame"){
                  dat[,numcol]<-discretize(data.frame(dat[,numcol]),method="quantile")
                }else{ print(tried)
                  stop("the numeric covariate could not be discretized, recode it into factor")
                }
              }else{
                tried<-try(discretize(dat[,numcol],method="quantile"),silent=TRUE)
                if(class(tried)=="data.frame"){
                  dat[,numcol]<-discretize(dat[,numcol],method="quantile")
                }else{ print(tried)
                  stop("at least one numeric covariate could not be discretized, recode it into factor")
                }
              }} 
            
          }
          #####Markov networks
          
          if(type=="mmpc"){
            #Common cause criterion (Waernbaum, de Luna and Richardson)
            ##Algorithm 1
            bmT<-matrix(c(rep("T",length(covarcol)),covars),ncol=2)
            blacklistT<-data.frame(bmT)
            names(blacklistT)<-c("from","to")
            lt1<-dat[,c(covarcol,Tcol)]
            res1<-mmpc(lt1,optimized=FALSE)
            ##Subset X.T
            XT<-res1$nodes$T$mb
            covarsT<-which(match(colnames(dat),XT)!="NA")
            bmQ<-matrix(c(rep("Y",length(covarsT)),XT),ncol=2)
            blacklistQ<-data.frame(bmQ)
            names(blacklistQ)<-c("from","to")
            lt2<-dat[which(dat[,Tcol]==1),c(covarsT,ycol)]
            ##Subset Q.1
            if(is.data.frame(lt2)==FALSE){
              Q1<-NULL
            }else{
              res2<-mmpc(lt2,optimized=FALSE)
              Q1<-res2$nodes$Y$mb
            }
            lt3<-dat[which(dat[,Tcol]==0),c(covarsT,ycol)]
            ##Subset Q.0
            if(is.data.frame(lt3)==FALSE){
              Q0<-NULL
            }else{
              res3<-mmpc(lt3,optimized=FALSE)
              Q0<-res3$nodes$Y$mb
            }
            ##Subset Q=Q.1UQ.0
            if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}
            
            
            ##Algortithm 2
            bmY<-matrix(c(rep("Y",length(covarcol)),covars),ncol=2)
            blacklistY<-data.frame(bmY)
            names(blacklistY)<-c("from","to")
            lt4<-dat[which(dat[,Tcol]==1),c(covarcol,ycol)]
            lt5<-dat[which(dat[,Tcol]==0),c(covarcol,ycol)]
            
            res4<-mmpc(lt4,optimized=FALSE)
            res5<-mmpc(lt5,optimized=FALSE)
            
            
            ##Subset X.1
            X1<-res4$nodes$Y$mb
            ##Subset X.0
            X0<-res5$nodes$Y$mb
            ##Subset X.Y
            if((length(X1)+length(X0)==0)){XY<-NULL}else{
              XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
            covars1<-which(match(colnames(dat),X1)!="NA")
            covars0<-which(match(colnames(dat),X0)!="NA")
            bm1<-matrix(c(rep("T",length(covars1)),X1),ncol=2)
            blacklist1<-data.frame(bm1)
            names(blacklist1)<-c("from","to")
            bm0<-matrix(c(rep("T",length(covars0)),X0),ncol=2)
            blacklist0<-data.frame(bm0)
            names(blacklist0)<-c("from","to")
            lt6<-dat[,c(covars1,Tcol)]
            lt7<-dat[,c(covars0,Tcol)]
            ##Subset Z.1
            if(is.data.frame(lt6)==FALSE){
              Z1<-NULL
            }else{
              res6<-mmpc(lt6,optimized=FALSE)
              Z1<-res6$nodes$T$mb
            }
            ##Subset Z.0
            if(is.data.frame(lt7)==FALSE){
              Z0<-NULL
            }else{
              res7<-mmpc(lt7,optimized=FALSE)
              Z0<-res7$nodes$T$mb
            }
            ##Subset Z=Z.1UZ.0
            if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}
            
            #Disjunctive cause criterion (VanderWeele and Shpitser)
            ##Subset X.D
            if((length(XY)+length(XT)==0)){XD<-NULL}else{
              XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
            }
            
            covarsXD<-which(match(colnames(dat),XD)!="NA")
            
            
          }
          
          
          #########Bayesian networks
          
          if(type=="mmhc"){
            #Common cause criterion (Waernbaum, de Luna and Richardson)
            ##Algorithm 1
            bmT<-matrix(c(rep("T",length(covarcol)),covars),ncol=2)
            blacklistT<-data.frame(bmT)
            names(blacklistT)<-c("from","to")
            lt1<-dat[,c(covarcol,Tcol)]
            res1<-mmhc(lt1,blacklist=blacklistT,optimized=FALSE, score="aic")
            ##Subset X.T
            XT<-res1$nodes$T$mb
            #print(XT)
            covarsT<-which(match(colnames(dat),XT)!="NA")
            bmQ<-matrix(c(rep("Y",length(covarsT)),XT),ncol=2)
            blacklistQ<-data.frame(bmQ)
            names(blacklistQ)<-c("from","to")
            lt2<-dat[which(dat[,Tcol]==1),c(covarsT,ycol)]
            ##Subset Q.1
            if(is.data.frame(lt2)==FALSE){
              Q1<-NULL
            }else{
              res2<-mmhc(lt2,blacklist=blacklistQ,optimized=FALSE, score="aic")
              Q1<-res2$nodes$Y$mb
            }
            lt3<-dat[which(dat[,Tcol]==0),c(covarsT,ycol)]
            ##Subset Q.0
            if(is.data.frame(lt3)==FALSE){
              Q0<-NULL
            }else{
              res3<-mmhc(lt3,blacklist=blacklistQ,optimized=FALSE, score="aic")
              Q0<-res3$nodes$Y$mb
            }
            ##Subset Q=Q.1UQ.0
            if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}
            
            
            ##Algortithm 2
            bmY<-matrix(c(rep("Y",length(covarcol)),covars),ncol=2)
            blacklistY<-data.frame(bmY)
            names(blacklistY)<-c("from","to")
            lt4<-dat[which(dat[,Tcol]==1),c(covarcol,ycol)]
            lt5<-dat[which(dat[,Tcol]==0),c(covarcol,ycol)]
            
            
            res4<-mmhc(lt4,blacklist=blacklistY,optimized=FALSE, score="aic")
            res5<-mmhc(lt5,blacklist=blacklistY,optimized=FALSE, score="aic")
            
            ##Subset X.1
            X1<-res4$nodes$Y$mb
            ##Subset X.0
            X0<-res5$nodes$Y$mb
            ##Subset X.Y
            if((length(X1)+length(X0)==0)){XY<-NULL}else{
              XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
            covars1<-which(match(colnames(dat),X1)!="NA")
            covars0<-which(match(colnames(dat),X0)!="NA")
            bm1<-matrix(c(rep("T",length(covars1)),X1),ncol=2)
            blacklist1<-data.frame(bm1)
            names(blacklist1)<-c("from","to")
            bm0<-matrix(c(rep("T",length(covars0)),X0),ncol=2)
            blacklist0<-data.frame(bm0)
            names(blacklist0)<-c("from","to")
            lt6<-dat[,c(covars1,Tcol)]
            lt7<-dat[,c(covars0,Tcol)]
            ##Subset Z.1
            if(is.data.frame(lt6)==FALSE){
              Z1<-NULL
            }else{
              res6<-mmhc(lt6,blacklist=blacklist1,optimized=FALSE, score="aic")
              Z1<-res6$nodes$T$mb
            }
            ##Subset Z.0
            if(is.data.frame(lt7)==FALSE){
              Z0<-NULL
            }else{
              res7<-mmhc(lt7,blacklist=blacklist0,optimized=FALSE, score="aic")
              Z0<-res7$nodes$T$mb
            }
            ##Subset Z=Z.1UZ.0
            if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}
            
            #Disjunctive cause criterion (VanderWeele and Shpitser)
            ##Subset X.D
            if((length(XY)+length(XT)==0)){XD<-NULL}else{
              XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
            }
            covarsXD<-which(match(colnames(dat),XD)!="NA")
            
            
          }
          
          ######LASSO
          
          if(type=="lasso"){
            #Common cause criterion (Waernbaum, de Luna and Richardson)
            ##Algorithm 1
            Y<-datbeta$T
            X<-datbeta[,covarcol]
            D<-dim(X)[2]
            
            Xo<-X
            if(length(numcol)>0){
              xnum<-X[,numcol]
              Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
              f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol], '^2)', sep=''), collapse=" + ")))
              f2<-as.formula(f1)
              x1<- model.matrix(f2, cbind(X,Y))[, -1]
              f3 <- as.formula(Y ~ .*.)
              x2<- model.matrix(f3, cbind(X,Y))[, -1]
              X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
              
              
            }else{
              Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
              
              f3 <- as.formula(Y ~ .*.)
              x2<- model.matrix(f3, data.frame(cbind(X,Y)))[, -1]
              X<-cbind(x2[,c(1:D)])
              
              
            }
            SL <-cov.sel.high.lasso(Y = Y, X = X)
            
            ##Subset X.T
            u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
            
            XT<-covars[is.na(match(covars,u))==FALSE]
            covarsT<-which(match(colnames(datbeta),XT)!="NA")
            
            lt2<-datbeta[which(datbeta[,Tcol]==1),c(covarsT,ycol)]
            
            ##Subset Q.1
            if(is.data.frame(lt2)==FALSE){
              Q1<-NULL
            }else{
              Y<-lt2$Y
              X<-lt2[,-dim(lt2)[2]]
              xnames<-names(lt2)[-dim(lt2)[2]]
              D<-length(xnames)
              if(D==1){Q1<-XT}else{
                Xo<-X
                numcol2<-which(lapply(X,class)=="numeric")
                
                if(length(numcol2)>0){
                  
                  xnum<-X[,numcol2]
                  Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                  
                  f1<-eval(paste("Y ~", paste(paste('I(', xnames[numcol2], '^2)', sep=''), collapse=" + ")))
                  f2<-as.formula(f1)
                  x1<- model.matrix(f2, cbind(X,Y))[, -1]
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
                  
                }else{
                  Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                  
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)])
                  
                  
                }
                
                
                SL <-cov.sel.high.lasso(Y = Y, X = X)
                
                
                u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
                Q1<-covars[is.na(match(covars,u))==FALSE]
                
                
              }
            }
            
            
            lt2<-datbeta[which(datbeta[,Tcol]==0),c(covarsT,ycol)]
            
            ##Subset Q.0
            if(is.data.frame(lt2)==FALSE){
              Q0<-NULL
            }else{
              Y<-lt2$Y
              X<-lt2[,-dim(lt2)[2]]
              xnames<-names(lt2)[-dim(lt2)[2]]
              D<-length(xnames)
              if(D==1){Q0<-XT}else{
                
                Xo<-X
                numcol2<-which(lapply(X,class)=="numeric")
                if(length(numcol2)>0){
                  xnum<-X[,numcol2]
                  Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                  f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol2], '^2)', sep=''), collapse=" + ")))
                  f2<-as.formula(f1)
                  x1<- model.matrix(f2, cbind(X,Y))[, -1]
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
                }else{
                  Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)])
                  
                  
                }
                
                
                
                SL <-cov.sel.high.lasso(Y = Y, X = X)
                
                
                u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
                Q0<-covars[is.na(match(covars,u))==FALSE]
                
              }
              
            }
            ##Subset Q=Q.1UQ.0
            if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}
            
            ##Algortithm 2
            lt4<-datbeta[which(datbeta[,Tcol]==1),c(covarcol,ycol)]
            Y<-lt4$Y
            X<-lt4[,covarcol]
            D<-dim(X)[2]
            
            Xo<-X
            if(length(numcol)>0){
              xnum<-X[,numcol]
              Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
              f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol], '^2)', sep=''), collapse=" + ")))
              f2<-as.formula(f1)
              x1<- model.matrix(f2, cbind(X,Y))[, -1]
              f3 <- as.formula(Y ~ .*.)
              x2<- model.matrix(f3, cbind(X,Y))[, -1]
              X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
            }else{
              Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
              f3 <- as.formula(Y ~ .*.)
              x2<- model.matrix(f3, cbind(X,Y))[, -1]
              X<-cbind(x2[,c(1:D)])
              
              
            }
            
            SL <-cov.sel.high.lasso(Y = Y, X = X)
            
            ##Subset X.1
            u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
            
            X1<-covars[is.na(match(covars,u))==FALSE]
            covarsX1<-which(match(colnames(datbeta),X1)!="NA")
            
            lt2<-datbeta[,c(covarsX1,Tcol)]
            
            ##Subset Z.1
            if(is.data.frame(lt2)==FALSE){
              Z1<-NULL
            }else{
              Y<-lt2$T
              X<-lt2[,-dim(lt2)[2]]
              xnames<-names(lt2)[-dim(lt2)[2]]
              D<-length(xnames)
              if(D==1){Z1<-X1}else{
                Xo<-X
                numcol2<-which(lapply(X,class)=="numeric")
                if(length(numcol2)>0){
                  xnum<-X[,numcol2]
                  Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                  f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol2], '^2)', sep=''), collapse=" + ")))
                  f2<-as.formula(f1)
                  x1<- model.matrix(f2, cbind(X,Y))[, -1]
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
                }else{
                  Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)])
                  
                  
                }
                
                SL <-cov.sel.high.lasso(Y = Y, X = X)
                
                
                u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
                Z1<-covars[is.na(match(covars,u))==FALSE]
                
              }
              
            }
            lt4<-datbeta[which(datbeta[,Tcol]==0),c(covarcol,ycol)]
            Y<-lt4$Y
            X<-lt4[,covarcol]
            D<-dim(X)[2]
            
            Xo<-X
            if(length(numcol)>0){
              xnum<-X[,numcol]
              Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
              
              f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol], '^2)', sep=''), collapse=" + ")))
              f2<-as.formula(f1)
              x1<- model.matrix(f2, cbind(X,Y))[, -1]
              f3 <- as.formula(Y ~ .*.)
              x2<- model.matrix(f3, cbind(X,Y))[, -1]
              X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
            }else{
              Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
              f3 <- as.formula(Y ~ .*.)
              x2<- model.matrix(f3, cbind(X,Y))[, -1]
              X<-cbind(x2[,c(1:D)])
              
              
            }
            
            SL <-cov.sel.high.lasso(Y = Y, X = X)
            
            ##Subset X.0
            u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
            X0<-covars[is.na(match(covars,u))==FALSE]
            covarsX0<-which(match(colnames(datbeta),X0)!="NA")
            
            lt2<-datbeta[,c(covarsX0,Tcol)]
            
            ##Subset Z.0
            if(is.data.frame(lt2)==FALSE){
              Z0<-NULL
            }else{
              Y<-lt2$T
              X<-lt2[,-dim(lt2)[2]]
              xnames<-names(lt2)[-dim(lt2)[2]]
              D<-length(xnames)
              if(D==1){Z0<-X0}else{
                Xo<-X
                numcol2<-which(lapply(X,class)=="numeric")
                if(length(numcol2)>0){
                  xnum<-X[,numcol2]
                  Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                  
                  f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol2], '^2)', sep=''), collapse=" + ")))
                  f2<-as.formula(f1)
                  x1<- model.matrix(f2, cbind(X,Y))[, -1]
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
                }else{
                  Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                  
                  f3 <- as.formula(Y ~ .*.)
                  x2<- model.matrix(f3, cbind(X,Y))[, -1]
                  X<-cbind(x2[,c(1:D)])
                  
                  
                }
                
                SL <-cov.sel.high.lasso(Y = Y, X = X)
                
                u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
                
                Z0<-covars[is.na(match(covars,u))==FALSE]
                
              }
            }
            
            ##Subset X.Y
            if((length(X1)+length(X0)==0)){XY<-NULL}else{
              XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
            
            
            ##Subset Z=Z.1UZ.0
            if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}
            
            #Disjunctive cause criterion (VanderWeele and Shpitser)
            ##Subset X.D
            if((length(XY)+length(XT)==0)){XD<-NULL}else{
              XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
            }
            covarsXD<-which(match(colnames(datbeta),XD)!="NA")
            
          }
          
          
          if(type=="rf"){
            #Common cause criterion (Waernbaum, de Luna and Richardson)
            ##Algorithm 1
            Y<-datbeta$T
            X<-datbeta[,covarcol]
            
            
            SL <-cov.sel.high.rf(Y = Y, X = X)
            
            ##Subset X.T
            XT<-colnames(X)[which(SL==TRUE)]
            covarsT<-which(match(colnames(datbeta),XT)!="NA")
            
            lt2<-datbeta[which(datbeta[,Tcol]==1),c(covarsT,ycol)]
            
            ##Subset Q.1
            if(is.data.frame(lt2)==FALSE){
              Q1<-NULL
            }else{
              Y<-lt2$Y
              X<-lt2[,-dim(lt2)[2]]
              if(is.data.frame(X)==FALSE){
                X<-data.frame(X)
                names(X)<-XT
              }
              SL <-cov.sel.high.rf(Y = Y, X = X)
              Q1<-colnames(X)[which(SL==TRUE)]
              
              
            }
            
            
            
            lt2<-datbeta[which(datbeta[,Tcol]==0),c(covarsT,ycol)]
            
            ##Subset Q.0
            if(is.data.frame(lt2)==FALSE){
              Q0<-NULL
            }else{
              Y<-lt2$Y
              X<-lt2[,-dim(lt2)[2]]
              if(is.data.frame(X)==FALSE){
                X<-data.frame(X)
                names(X)<-XT
              }
              SL <-cov.sel.high.rf(Y = Y, X = X)
              Q0<-colnames(X)[which(SL==TRUE)]
              
            }
            
            
            ##Subset Q=Q.1UQ.0
            if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}
            
            ##Algortithm 2
            lt4<-datbeta[which(datbeta[,Tcol]==1),c(covarcol,ycol)]
            Y<-lt4$Y
            X<-lt4[,covarcol]
            
            SL <-cov.sel.high.rf(Y = Y, X = X)
            
            
            ##Subset X.1
            
            X1<-colnames(X)[which(SL==TRUE)]
            covarsX1<-which(match(colnames(datbeta),X1)!="NA")
            
            lt2<-datbeta[,c(covarsX1,Tcol)]
            
            ##Subset Z.1
            if(is.data.frame(lt2)==FALSE){
              Z1<-NULL
            }else{
              Y<-lt2$T
              X<-lt2[,-dim(lt2)[2]]
              
              if(is.data.frame(X)==FALSE){
                X<-data.frame(X)
                names(X)<-X1
              }
              SL <-cov.sel.high.rf(Y = Y, X = X)
              Z1<-colnames(X)[which(SL==TRUE)]
              
            }
            
            
            lt4<-datbeta[which(datbeta[,Tcol]==0),c(covarcol,ycol)]
            Y<-lt4$Y
            X<-lt4[,covarcol]
            SL <-cov.sel.high.rf(Y = Y, X = X)
            
            ##Subset X.0
            X0<-colnames(X)[which(SL==TRUE)]
            covarsX0<-which(match(colnames(datbeta),X0)!="NA")
            
            lt2<-datbeta[,c(covarsX0,Tcol)]
            
            ##Subset Z.0
            if(is.data.frame(lt2)==FALSE){
              Z0<-NULL
            }else{
              Y<-lt2$T
              X<-lt2[,-dim(lt2)[2]]
              if(is.data.frame(X)==FALSE){
                X<-data.frame(X)
                names(X)<-X0
              }
              SL <-cov.sel.high.rf(Y = Y, X = X)
              Z0<-colnames(X)[which(SL==TRUE)]
              
            }
            
            
            ##Subset X.Y
            if((length(X1)+length(X0)==0)){XY<-NULL}else{
              XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
            
            
            ##Subset Z=Z.1UZ.0
            if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}
            
            #Disjunctive cause criterion (VanderWeele and Shpitser)
            ##Subset X.D
            if((length(XY)+length(XT)==0)){XD<-NULL}else{
              XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
            }
            covarsXD<-which(match(colnames(datbeta),XD)!="NA")
            
          }
          
          #Subset cardinalities
          cardXT<-length(XT)
          cardQ<-length(Q)
          cardXY<-length(XY)
          cardZ<-length(Z)
          cardXD<-length(XD)
          cards<-data.frame(X.T=cardXT,Q=cardQ,X.Y=cardXY,Z=cardZ,X.D=cardXD)
          
          #ATE estimate via propensity score matching
          if(betahat==TRUE){
            datbeta[,Tcol]<-as.numeric(datbeta[,Tcol])-1
            #Pre-treatment criterion
            f1<-as.formula(paste("T~", paste(paste(colnames(datbeta)[covarcol]), collapse= "+")))
            ps<-glm(f1,family=binomial,data=datbeta)$fitted
            betahatX<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1, replace=TRUE)
            betahatXest<-betahatX$est
            betahatXse<-betahatX$se
            
            #Common cause criterion
            ##Algorithm 1
            ##Subset X.T
            if(length(XT)==0){
              betahatXTest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
              betahatXTse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
              
            }else{
              f1<-as.formula(paste("T~", paste(paste(XT), collapse= "+")))
              ps<-glm(f1,family=binomial,data=datbeta)$fitted
              betahatXT<- Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
              betahatXTest<-betahatXT$est
              betahatXTse<-betahatXT$se
              
            }
            ##Subset Q
            if(length(Q)==0){   betahatQest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
            betahatQse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
            
            }else{
              f1<-as.formula(paste("T~", paste(paste(Q), collapse= "+")))
              ps<-glm(f1,family=binomial,data=datbeta)$fitted
              betahatQ<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
              betahatQest<-betahatQ$est
              betahatQse<-betahatQ$se
              
            }
            
            #Algorithm 2
            ##Subset X.Y
            if(length(XY)==0){   betahatXYest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
            betahatXYse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
            }else{
              f1<-as.formula(paste("T~", paste(paste(XY), collapse= "+")))
              ps<-glm(f1,family=binomial,data=datbeta)$fitted
              betahatXY<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
              betahatXYest<-betahatXY$est
              betahatXYse<-betahatXY$se
            }
            
            ##Subset Z
            if(length(Z)==0){   betahatZest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
            betahatZse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
            }else{
              f1<-as.formula(paste("T~", paste(paste(Z), collapse= "+")))
              ps<-glm(f1,family=binomial,data=datbeta)$fitted
              betahatZ<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,   replace=TRUE)
              betahatZest<-betahatZ$est
              betahatZse<-betahatZ$se
              
            }
            
            #Disjunctive cause criterion
            ##Subset X.D
            if(length(XD)==0){   betahatXDest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
            betahatXDse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
            }else{
              f1<-as.formula(paste("T~", paste(paste(XD), collapse= "+")))
              ps<-glm(f1,family=binomial,data=datbeta)$fitted
              betahatXD<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
              betahatXDest<-betahatXD$est
              betahatXDse<-betahatXD$se
              
            }
            
            
            
            
            betahats<-data.frame(X=betahatXest,X.T=betahatXTest,Q=betahatQest,X.Y=betahatXYest,Z=betahatZest,X.D=betahatXDest)
            betahatsse<-data.frame(X=betahatXse,X.T=betahatXTse,Q=betahatQse,X.Y=betahatXYse,Z=betahatZse,X.D=betahatXDse)
          }else{betahats<-betahatsse<-NULL}
          
      list(X.T=XT,Q.0=Q0,Q.1=Q1,Q=Q,X.0=X0,X.1=X1,X.Y=XY,Z.0=Z0,Z.1=Z1,Z=Z,X.TY=XD,cardinalities=cards, est=betahats, se=betahatsse, N=N, Setting=Setting, rep=Rep, Models=Models,type=type, varnames=varnames)

      }
    }else{

      for(i in 1:Rep){
        if(Simulate==TRUE){
          dat<-cov.sel.high.sim(N, Setting, Rep, Models)$dat
          X<-dat[,-c((dim(dat)[2])-1,dim(dat)[2])]
          Y<-dat[,((dim(dat)[2])-1)]
          T<-dat[,dim(dat)[2]]
        }
        varnames<-colnames(dat)
        
        covarcol<-1:(dim(X)[2])
        ycol<-(dim(X)[2])+1
        Tcol<-ycol+1
        if(class(dat[,Tcol])!="factor"){
          dat[,Tcol]<-factor(dat[,Tcol])
        }
        
        datbeta<-dat
        if(length(unique(dat[,ycol]))==2){
          datbeta[,ycol]<-factor(datbeta[,ycol])
        }
        
        
        
        
        
        
        numcol<-which(lapply(dat[covarcol],class)=="numeric")
        
        
        covars<-colnames(dat[,covarcol])
        
        
        if(type=="mmpc" || type=="mmhc"){
          if(length(unique(dat[,ycol]))>2){
            tried<-try(discretize(data.frame(dat[,ycol]),method="quantile"),silent=TRUE)
            if(class(tried)=="data.frame"){
              dat[,ycol]<-discretize(data.frame(dat[,ycol]),method="quantile")
            }else{
              print(tried)
              stop("the numeric outcome could not be discretized, recode it into factor")
              
            }
          }else{dat[,ycol]<-factor(dat[,ycol])
          }
          
          if(length(numcol)>0){
            if(length(numcol)==1){
              tried<-try(discretize(data.frame(dat[,numcol]),method="quantile"),silent=TRUE)
              if(class(tried)=="data.frame"){
                dat[,numcol]<-discretize(data.frame(dat[,numcol]),method="quantile")
              }else{ print(tried)
                stop("the numeric covariate could not be discretized, recode it into factor")
              }
            }else{
              tried<-try(discretize(dat[,numcol],method="quantile"),silent=TRUE)
              if(class(tried)=="data.frame"){
                dat[,numcol]<-discretize(dat[,numcol],method="quantile")
              }else{ print(tried)
                stop("at least one numeric covariate could not be discretized, recode it into factor")
              }
            }} 
          
        }
        #####Markov networks

        if(type=="mmpc"){
          #Common cause criterion (Waernbaum, de Luna and Richardson)
          ##Algorithm 1
          bmT<-matrix(c(rep("T",length(covarcol)),covars),ncol=2)
          blacklistT<-data.frame(bmT)
          names(blacklistT)<-c("from","to")
          lt1<-dat[,c(covarcol,Tcol)]
          res1<-mmpc(lt1,optimized=FALSE)
          ##Subset X.T
          XT<-res1$nodes$T$mb
          covarsT<-which(match(colnames(dat),XT)!="NA")
          bmQ<-matrix(c(rep("Y",length(covarsT)),XT),ncol=2)
          blacklistQ<-data.frame(bmQ)
          names(blacklistQ)<-c("from","to")
          lt2<-dat[which(dat[,Tcol]==1),c(covarsT,ycol)]
          ##Subset Q.1
          if(is.data.frame(lt2)==FALSE){
            Q1<-NULL
          }else{
            res2<-mmpc(lt2,optimized=FALSE)
            Q1<-res2$nodes$Y$mb
          }
          lt3<-dat[which(dat[,Tcol]==0),c(covarsT,ycol)]
          ##Subset Q.0
          if(is.data.frame(lt3)==FALSE){
            Q0<-NULL
          }else{
            res3<-mmpc(lt3,optimized=FALSE)
            Q0<-res3$nodes$Y$mb
          }
          ##Subset Q=Q.1UQ.0
          if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}


          ##Algortithm 2
          bmY<-matrix(c(rep("Y",length(covarcol)),covars),ncol=2)
          blacklistY<-data.frame(bmY)
          names(blacklistY)<-c("from","to")
          lt4<-dat[which(dat[,Tcol]==1),c(covarcol,ycol)]
          lt5<-dat[which(dat[,Tcol]==0),c(covarcol,ycol)]

          res4<-mmpc(lt4,optimized=FALSE)
          res5<-mmpc(lt5,optimized=FALSE)


          ##Subset X.1
          X1<-res4$nodes$Y$mb
          ##Subset X.0
          X0<-res5$nodes$Y$mb
          ##Subset X.Y
          if((length(X1)+length(X0)==0)){XY<-NULL}else{
          XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
          covars1<-which(match(colnames(dat),X1)!="NA")
          covars0<-which(match(colnames(dat),X0)!="NA")
          bm1<-matrix(c(rep("T",length(covars1)),X1),ncol=2)
          blacklist1<-data.frame(bm1)
          names(blacklist1)<-c("from","to")
          bm0<-matrix(c(rep("T",length(covars0)),X0),ncol=2)
          blacklist0<-data.frame(bm0)
          names(blacklist0)<-c("from","to")
          lt6<-dat[,c(covars1,Tcol)]
          lt7<-dat[,c(covars0,Tcol)]
          ##Subset Z.1
          if(is.data.frame(lt6)==FALSE){
            Z1<-NULL
          }else{
            res6<-mmpc(lt6,optimized=FALSE)
            Z1<-res6$nodes$T$mb
          }
          ##Subset Z.0
          if(is.data.frame(lt7)==FALSE){
            Z0<-NULL
          }else{
            res7<-mmpc(lt7,optimized=FALSE)
            Z0<-res7$nodes$T$mb
          }
          ##Subset Z=Z.1UZ.0
          if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}

          #Disjunctive cause criterion (VanderWeele and Shpitser)
          ##Subset X.D
          if((length(XY)+length(XT)==0)){XD<-NULL}else{
            XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
          }
          
          covarsXD<-which(match(colnames(dat),XD)!="NA")


        }


        #########Bayesian networks

        if(type=="mmhc"){
          #Common cause criterion (Waernbaum, de Luna and Richardson)
          ##Algorithm 1
          bmT<-matrix(c(rep("T",length(covarcol)),covars),ncol=2)
          blacklistT<-data.frame(bmT)
          names(blacklistT)<-c("from","to")
          lt1<-dat[,c(covarcol,Tcol)]
          res1<-mmhc(lt1,blacklist=blacklistT,optimized=FALSE, score="aic")
          ##Subset X.T
          XT<-res1$nodes$T$mb
          #print(XT)
          covarsT<-which(match(colnames(dat),XT)!="NA")
          bmQ<-matrix(c(rep("Y",length(covarsT)),XT),ncol=2)
          blacklistQ<-data.frame(bmQ)
          names(blacklistQ)<-c("from","to")
          lt2<-dat[which(dat[,Tcol]==1),c(covarsT,ycol)]
          ##Subset Q.1
          if(is.data.frame(lt2)==FALSE){
            Q1<-NULL
          }else{
            res2<-mmhc(lt2,blacklist=blacklistQ,optimized=FALSE, score="aic")
            Q1<-res2$nodes$Y$mb
          }
          lt3<-dat[which(dat[,Tcol]==0),c(covarsT,ycol)]
          ##Subset Q.0
          if(is.data.frame(lt3)==FALSE){
            Q0<-NULL
          }else{
            res3<-mmhc(lt3,blacklist=blacklistQ,optimized=FALSE, score="aic")
            Q0<-res3$nodes$Y$mb
          }
          ##Subset Q=Q.1UQ.0
          if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}


          ##Algortithm 2
          bmY<-matrix(c(rep("Y",length(covarcol)),covars),ncol=2)
          blacklistY<-data.frame(bmY)
          names(blacklistY)<-c("from","to")
          lt4<-dat[which(dat[,Tcol]==1),c(covarcol,ycol)]
          lt5<-dat[which(dat[,Tcol]==0),c(covarcol,ycol)]


          res4<-mmhc(lt4,blacklist=blacklistY,optimized=FALSE, score="aic")
          res5<-mmhc(lt5,blacklist=blacklistY,optimized=FALSE, score="aic")

          ##Subset X.1
          X1<-res4$nodes$Y$mb
          ##Subset X.0
          X0<-res5$nodes$Y$mb
          ##Subset X.Y
          if((length(X1)+length(X0)==0)){XY<-NULL}else{
            XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
          covars1<-which(match(colnames(dat),X1)!="NA")
          covars0<-which(match(colnames(dat),X0)!="NA")
          bm1<-matrix(c(rep("T",length(covars1)),X1),ncol=2)
          blacklist1<-data.frame(bm1)
          names(blacklist1)<-c("from","to")
          bm0<-matrix(c(rep("T",length(covars0)),X0),ncol=2)
          blacklist0<-data.frame(bm0)
          names(blacklist0)<-c("from","to")
          lt6<-dat[,c(covars1,Tcol)]
          lt7<-dat[,c(covars0,Tcol)]
          ##Subset Z.1
          if(is.data.frame(lt6)==FALSE){
            Z1<-NULL
          }else{
            res6<-mmhc(lt6,blacklist=blacklist1,optimized=FALSE, score="aic")
            Z1<-res6$nodes$T$mb
          }
          ##Subset Z.0
          if(is.data.frame(lt7)==FALSE){
            Z0<-NULL
          }else{
            res7<-mmhc(lt7,blacklist=blacklist0,optimized=FALSE, score="aic")
            Z0<-res7$nodes$T$mb
          }
          ##Subset Z=Z.1UZ.0
          if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}

          #Disjunctive cause criterion (VanderWeele and Shpitser)
          ##Subset X.D
          if((length(XY)+length(XT)==0)){XD<-NULL}else{
            XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
          }
          covarsXD<-which(match(colnames(dat),XD)!="NA")


        }

        ######LASSO

        if(type=="lasso"){
          #Common cause criterion (Waernbaum, de Luna and Richardson)
          ##Algorithm 1
          Y<-datbeta$T
          X<-datbeta[,covarcol]
          D<-dim(X)[2]
          
          Xo<-X
          if(length(numcol)>0){
            xnum<-X[,numcol]
            Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
            f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol], '^2)', sep=''), collapse=" + ")))
            f2<-as.formula(f1)
            x1<- model.matrix(f2, cbind(X,Y))[, -1]
            f3 <- as.formula(Y ~ .*.)
            x2<- model.matrix(f3, cbind(X,Y))[, -1]
            X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
            
            
          }else{
            Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
            
            f3 <- as.formula(Y ~ .*.)
            x2<- model.matrix(f3, data.frame(cbind(X,Y)))[, -1]
            X<-cbind(x2[,c(1:D)])
            
            
          }
          SL <-cov.sel.high.lasso(Y = Y, X = X)
          
          ##Subset X.T
          u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
          
          XT<-covars[is.na(match(covars,u))==FALSE]
          covarsT<-which(match(colnames(datbeta),XT)!="NA")
          
          lt2<-datbeta[which(datbeta[,Tcol]==1),c(covarsT,ycol)]
          
          ##Subset Q.1
          if(is.data.frame(lt2)==FALSE){
            Q1<-NULL
          }else{
            Y<-lt2$Y
            X<-lt2[,-dim(lt2)[2]]
            xnames<-names(lt2)[-dim(lt2)[2]]
            D<-length(xnames)
            if(D==1){Q1<-XT}else{
              Xo<-X
              numcol2<-which(lapply(X,class)=="numeric")
              
              if(length(numcol2)>0){
                
                xnum<-X[,numcol2]
                Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                
                f1<-eval(paste("Y ~", paste(paste('I(', xnames[numcol2], '^2)', sep=''), collapse=" + ")))
                f2<-as.formula(f1)
                x1<- model.matrix(f2, cbind(X,Y))[, -1]
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
                
              }else{
                Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)])
                
                
              }
              
              
              SL <-cov.sel.high.lasso(Y = Y, X = X)
              
              
              u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
              Q1<-covars[is.na(match(covars,u))==FALSE]
              
              
            }
          }
          
          
          lt2<-datbeta[which(datbeta[,Tcol]==0),c(covarsT,ycol)]
          
          ##Subset Q.0
          if(is.data.frame(lt2)==FALSE){
            Q0<-NULL
          }else{
            Y<-lt2$Y
            X<-lt2[,-dim(lt2)[2]]
            xnames<-names(lt2)[-dim(lt2)[2]]
            D<-length(xnames)
            if(D==1){Q0<-XT}else{
              
              Xo<-X
              numcol2<-which(lapply(X,class)=="numeric")
              if(length(numcol2)>0){
                xnum<-X[,numcol2]
                Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol2], '^2)', sep=''), collapse=" + ")))
                f2<-as.formula(f1)
                x1<- model.matrix(f2, cbind(X,Y))[, -1]
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
              }else{
                Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)])
                
                
              }
              
              
              
              SL <-cov.sel.high.lasso(Y = Y, X = X)
              
              
              u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
              Q0<-covars[is.na(match(covars,u))==FALSE]
              
            }
            
          }
          ##Subset Q=Q.1UQ.0
          if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}
          
          ##Algortithm 2
          lt4<-datbeta[which(datbeta[,Tcol]==1),c(covarcol,ycol)]
          Y<-lt4$Y
          X<-lt4[,covarcol]
          D<-dim(X)[2]
          
          Xo<-X
          if(length(numcol)>0){
            xnum<-X[,numcol]
            Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
            f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol], '^2)', sep=''), collapse=" + ")))
            f2<-as.formula(f1)
            x1<- model.matrix(f2, cbind(X,Y))[, -1]
            f3 <- as.formula(Y ~ .*.)
            x2<- model.matrix(f3, cbind(X,Y))[, -1]
            X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
          }else{
            Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
            f3 <- as.formula(Y ~ .*.)
            x2<- model.matrix(f3, cbind(X,Y))[, -1]
            X<-cbind(x2[,c(1:D)])
            
            
          }
          
          SL <-cov.sel.high.lasso(Y = Y, X = X)
          
          ##Subset X.1
          u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
          
          X1<-covars[is.na(match(covars,u))==FALSE]
          covarsX1<-which(match(colnames(datbeta),X1)!="NA")
          
          lt2<-datbeta[,c(covarsX1,Tcol)]
          
          ##Subset Z.1
          if(is.data.frame(lt2)==FALSE){
            Z1<-NULL
          }else{
            Y<-lt2$T
            X<-lt2[,-dim(lt2)[2]]
            xnames<-names(lt2)[-dim(lt2)[2]]
            D<-length(xnames)
            if(D==1){Z1<-X1}else{
              Xo<-X
              numcol2<-which(lapply(X,class)=="numeric")
              if(length(numcol2)>0){
                xnum<-X[,numcol2]
                Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol2], '^2)', sep=''), collapse=" + ")))
                f2<-as.formula(f1)
                x1<- model.matrix(f2, cbind(X,Y))[, -1]
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
              }else{
                Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)])
                
                
              }
              
              SL <-cov.sel.high.lasso(Y = Y, X = X)
              
              
              u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
              Z1<-covars[is.na(match(covars,u))==FALSE]
              
            }
            
          }
          lt4<-datbeta[which(datbeta[,Tcol]==0),c(covarcol,ycol)]
          Y<-lt4$Y
          X<-lt4[,covarcol]
          D<-dim(X)[2]
          
          Xo<-X
          if(length(numcol)>0){
            xnum<-X[,numcol]
            Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
            
            f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol], '^2)', sep=''), collapse=" + ")))
            f2<-as.formula(f1)
            x1<- model.matrix(f2, cbind(X,Y))[, -1]
            f3 <- as.formula(Y ~ .*.)
            x2<- model.matrix(f3, cbind(X,Y))[, -1]
            X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
          }else{
            Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
            f3 <- as.formula(Y ~ .*.)
            x2<- model.matrix(f3, cbind(X,Y))[, -1]
            X<-cbind(x2[,c(1:D)])
            
            
          }
          
          SL <-cov.sel.high.lasso(Y = Y, X = X)
          
          ##Subset X.0
          u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
          X0<-covars[is.na(match(covars,u))==FALSE]
          covarsX0<-which(match(colnames(datbeta),X0)!="NA")
          
          lt2<-datbeta[,c(covarsX0,Tcol)]
          
          ##Subset Z.0
          if(is.data.frame(lt2)==FALSE){
            Z0<-NULL
          }else{
            Y<-lt2$T
            X<-lt2[,-dim(lt2)[2]]
            xnames<-names(lt2)[-dim(lt2)[2]]
            D<-length(xnames)
            if(D==1){Z0<-X0}else{
              Xo<-X
              numcol2<-which(lapply(X,class)=="numeric")
              if(length(numcol2)>0){
                xnum<-X[,numcol2]
                Xvars<-t(rbind(matrix(c(names(X),names(xnum),rep(0,length(c(names(X),names(xnum))))),ncol=2),  t(combn(names(X),2))))
                
                f1<-eval(paste("Y ~", paste(paste('I(', names(X)[numcol2], '^2)', sep=''), collapse=" + ")))
                f2<-as.formula(f1)
                x1<- model.matrix(f2, cbind(X,Y))[, -1]
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)],x1,x2[,-c(1:D)])
              }else{
                Xvars<-t(rbind(matrix(c(names(X),rep(0,length(names(X)))),ncol=2),  t(combn(names(X),2))))
                
                f3 <- as.formula(Y ~ .*.)
                x2<- model.matrix(f3, cbind(X,Y))[, -1]
                X<-cbind(x2[,c(1:D)])
                
                
              }
              
              SL <-cov.sel.high.lasso(Y = Y, X = X)
              
              u<-unique(c(unique(Xvars[1,which(SL==TRUE)]),unique(Xvars[2,which(SL==TRUE)])))
              
              Z0<-covars[is.na(match(covars,u))==FALSE]
              
            }
          }
          
          ##Subset X.Y
          if((length(X1)+length(X0)==0)){XY<-NULL}else{
            XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
          
          
          ##Subset Z=Z.1UZ.0
          if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}
          
          #Disjunctive cause criterion (VanderWeele and Shpitser)
          ##Subset X.D
          if((length(XY)+length(XT)==0)){XD<-NULL}else{
            XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
          }
          covarsXD<-which(match(colnames(datbeta),XD)!="NA")
          
        }
        

        if(type=="rf"){
          #Common cause criterion (Waernbaum, de Luna and Richardson)
          ##Algorithm 1
          Y<-datbeta$T
          X<-datbeta[,covarcol]


          SL <-cov.sel.high.rf(Y = Y, X = X)

          ##Subset X.T
          XT<-colnames(X)[which(SL==TRUE)]
          covarsT<-which(match(colnames(datbeta),XT)!="NA")

          lt2<-datbeta[which(datbeta[,Tcol]==1),c(covarsT,ycol)]

          ##Subset Q.1
          if(is.data.frame(lt2)==FALSE){
            Q1<-NULL
          }else{
            Y<-lt2$Y
            X<-lt2[,-dim(lt2)[2]]
            if(is.data.frame(X)==FALSE){
              X<-data.frame(X)
              names(X)<-XT
            }
            SL <-cov.sel.high.rf(Y = Y, X = X)
            Q1<-colnames(X)[which(SL==TRUE)]


          }



          lt2<-datbeta[which(datbeta[,Tcol]==0),c(covarsT,ycol)]

          ##Subset Q.0
          if(is.data.frame(lt2)==FALSE){
            Q0<-NULL
          }else{
            Y<-lt2$Y
            X<-lt2[,-dim(lt2)[2]]
            if(is.data.frame(X)==FALSE){
              X<-data.frame(X)
              names(X)<-XT
            }
            SL <-cov.sel.high.rf(Y = Y, X = X)
            Q0<-colnames(X)[which(SL==TRUE)]

          }


          ##Subset Q=Q.1UQ.0
          if((length(Q1)+length(Q0)==0)){Q<-NULL}else{Q<-unique(c(Q1,Q0))[order(unique(c(Q1,Q0)))]}

          ##Algortithm 2
          lt4<-datbeta[which(datbeta[,Tcol]==1),c(covarcol,ycol)]
          Y<-lt4$Y
          X<-lt4[,covarcol]

          SL <-cov.sel.high.rf(Y = Y, X = X)


          ##Subset X.1

          X1<-colnames(X)[which(SL==TRUE)]
          covarsX1<-which(match(colnames(datbeta),X1)!="NA")

          lt2<-datbeta[,c(covarsX1,Tcol)]

          ##Subset Z.1
          if(is.data.frame(lt2)==FALSE){
            Z1<-NULL
          }else{
            Y<-lt2$T
            X<-lt2[,-dim(lt2)[2]]

            if(is.data.frame(X)==FALSE){
              X<-data.frame(X)
              names(X)<-X1
            }
            SL <-cov.sel.high.rf(Y = Y, X = X)
            Z1<-colnames(X)[which(SL==TRUE)]

          }


          lt4<-datbeta[which(datbeta[,Tcol]==0),c(covarcol,ycol)]
          Y<-lt4$Y
          X<-lt4[,covarcol]
          SL <-cov.sel.high.rf(Y = Y, X = X)

          ##Subset X.0
          X0<-colnames(X)[which(SL==TRUE)]
          covarsX0<-which(match(colnames(datbeta),X0)!="NA")

          lt2<-datbeta[,c(covarsX0,Tcol)]

          ##Subset Z.0
          if(is.data.frame(lt2)==FALSE){
            Z0<-NULL
          }else{
            Y<-lt2$T
            X<-lt2[,-dim(lt2)[2]]
            if(is.data.frame(X)==FALSE){
              X<-data.frame(X)
              names(X)<-X0
            }
            SL <-cov.sel.high.rf(Y = Y, X = X)
            Z0<-colnames(X)[which(SL==TRUE)]

          }


          ##Subset X.Y
          if((length(X1)+length(X0)==0)){XY<-NULL}else{
            XY<-unique(c(X1,X0))[order(unique(c(X1,X0)))]}
          

          ##Subset Z=Z.1UZ.0
          if((length(Z1)+length(Z0)==0)){Z<-NULL}else{Z<-unique(c(Z1,Z0))[order(unique(c(Z1,Z0)))]}

          #Disjunctive cause criterion (VanderWeele and Shpitser)
          ##Subset X.D
          if((length(XY)+length(XT)==0)){XD<-NULL}else{
            XD<-unique(c(XY,XT))[order(unique(c(XY,XT)))]
          }
          covarsXD<-which(match(colnames(datbeta),XD)!="NA")

        }

        #Subset cardinalities
        cardXT<-length(XT)
        cardQ<-length(Q)
        cardXY<-length(XY)
        cardZ<-length(Z)
        cardXD<-length(XD)
        cards<-data.frame(X.T=cardXT,Q=cardQ,X.Y=cardXY,Z=cardZ,X.D=cardXD)

        #ATE estimate via propensity score matching
        if(betahat==TRUE){
          datbeta[,Tcol]<-as.numeric(datbeta[,Tcol])-1
          #Pre-treatment criterion
          f1<-as.formula(paste("T~", paste(paste(colnames(datbeta)[covarcol]), collapse= "+")))
          ps<-glm(f1,family=binomial,data=datbeta)$fitted
          betahatX<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1, replace=TRUE)
          betahatXest<-betahatX$est
          betahatXse<-betahatX$se

          #Common cause criterion
          ##Algorithm 1
          ##Subset X.T
          if(length(XT)==0){
            betahatXTest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
            betahatXTse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))

          }else{
            f1<-as.formula(paste("T~", paste(paste(XT), collapse= "+")))
            ps<-glm(f1,family=binomial,data=datbeta)$fitted
            betahatXT<- Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
            betahatXTest<-betahatXT$est
            betahatXTse<-betahatXT$se

          }
          ##Subset Q
          if(length(Q)==0){   betahatQest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
          betahatQse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))

          }else{
            f1<-as.formula(paste("T~", paste(paste(Q), collapse= "+")))
            ps<-glm(f1,family=binomial,data=datbeta)$fitted
            betahatQ<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
            betahatQest<-betahatQ$est
            betahatQse<-betahatQ$se

          }

          #Algorithm 2
          ##Subset X.Y
          if(length(XY)==0){   betahatXYest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
          betahatXYse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
          }else{
            f1<-as.formula(paste("T~", paste(paste(XY), collapse= "+")))
            ps<-glm(f1,family=binomial,data=datbeta)$fitted
            betahatXY<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
            betahatXYest<-betahatXY$est
            betahatXYse<-betahatXY$se
          }

          ##Subset Z
          if(length(Z)==0){   betahatZest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
          betahatZse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
          }else{
            f1<-as.formula(paste("T~", paste(paste(Z), collapse= "+")))
            ps<-glm(f1,family=binomial,data=datbeta)$fitted
            betahatZ<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,   replace=TRUE)
            betahatZest<-betahatZ$est
            betahatZse<-betahatZ$se

          }

          #Disjunctive cause criterion
          ##Subset X.D
          if(length(XD)==0){   betahatXDest<-mean(datbeta[which(datbeta[,Tcol]==1),ycol])-mean(datbeta[which(datbeta[,Tcol]==0),ycol])
          betahatXDse<-sqrt((var(datbeta[which(datbeta[,Tcol]==1),ycol])/length(which(datbeta[,Tcol]==1))+var(datbeta[which(datbeta[,Tcol]==0),ycol]))/length(which(datbeta[,Tcol]==0)))
          }else{
            f1<-as.formula(paste("T~", paste(paste(XD), collapse= "+")))
            ps<-glm(f1,family=binomial,data=datbeta)$fitted
            betahatXD<-Match(Y=datbeta[,ycol], Tr=datbeta[,Tcol], X=as.matrix(ps,ncol=1), estimand = "ATE", M = 1,replace=TRUE)
            betahatXDest<-betahatXD$est
            betahatXDse<-betahatXD$se

          }




          betahats<-data.frame(X=betahatXest,X.T=betahatXTest,Q=betahatQest,X.Y=betahatXYest,Z=betahatZest,X.D=betahatXDest)
          betahatsse<-data.frame(X=betahatXse,X.T=betahatXTse,Q=betahatQse,X.Y=betahatXYse,Z=betahatZse,X.D=betahatXDse)
        }else{betahats<-betahatsse<-NULL}

        reslist[[i]]<-list(X.T=XT,Q.0=Q0,Q.1=Q1,Q=Q,X.0=X0,X.1=X1,X.Y=XY,Z.0=Z0,Z.1=Z1,Z=Z,X.TY=XD,cardinalities=cards, est=betahats, se=betahatsse, N=N, Setting=Setting, rep=Rep, Models=Models,type=type, varnames=varnames)
      }


    }

    l<-list(reslist=reslist)
    if(Simulate==TRUE){invisible(return(l[[1]]))}else{invisible(return(l[[1]][[1]]))}
}


