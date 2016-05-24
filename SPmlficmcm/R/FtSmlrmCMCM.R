FtSmlrmCMCM <-
function(fl,N,theta,beta,interc,vpo,vprob,vcorr){
                # Arguments function
                # fl : la formule 
                # N : l echantillon total 
                # theta : la probabilite de lallele mineur du genotype de de la mere 
                # beta : les parametre du modele de regression logistic sans intercepte 
                # prev est la prevalence de l evenement
                # vpo: is vector to indicated the position of mother genotype and chil genotype
                varz0<-all.vars(fl)
                varz1<-varz0[-1];
                varz2<-varz1[-vpo];n<-length(varz2)
                if(length(vprob)!=n){print(vprob)
                                      stop("length of probabilities vector is not equal to number of factors")} 
                if(max(vprob)>1){print(vprob)
                                      stop("the probabilities superior 1")}
                if(min(vprob)<0){print(vprob)
                                      stop("The probabilities inferior 0")} 
                if(length(vcorr)!=n){print(vcorr)
                                    stop("length of dependency coefficient vector is not equal to number of factors")
                                    }else{
                obs<-c(1:N)
                vecp=c((1-theta)**2,2*theta*(1-theta),theta**2)
                vecp0=c(1-theta,theta,0)
                vecp1=c((1-theta)/2,1/2,theta/2)
                vecp2=c(0,1-theta,theta)
                # genotype de la mere 
                gm<-sample(c(0,1,2),N,replace=TRUE,prob=vecp)
                # la fonction qui genere le genotype de enfant sachant celui de la mere 
                fgc1<-Vectorize(function(m){
                     if(m==0){gc<-sample(c(0,1,2),1,replace=TRUE,prob=vecp0)}
                     if(m==1){gc<-sample(c(0,1,2),1,replace=TRUE,prob=vecp1)}
                     if(m==2){gc<-sample(c(0,1,2),1,replace=TRUE,prob=vecp2)}
                     return(gc)},"m")
                # child genotype
                gc<-sapply(gm,fgc1)
                # the equation terms
                # la variable de l environement 
                dat<-NULL;
                for(i in 1:n){
                              dat<-cbind(dat,rbinom(length(gm),c(0,1),prob=vprob[i])+vcorr[i]*gm)
                              }
                # nom des varz
                nva<-c(varz0[1],varz1[-vpo],varz1[vpo[1]],varz1[vpo[2]])
                Y2<-rep(1,length(gm))
                # la base de donnees
                dat1<-data.frame(Y2,dat,gm,gc);
                names(dat1)<-nva;
                vx <- model.matrix(fl,data=dat1)
                # la fonction logistique
                logistic <- function(x){
                            return(1/(1+exp(-x)))
                            }
                # intercepte
                beta1=c(interc,beta)
                # la fonction qui l inverse du logit
                thet1 <- function(x){ 
                                    logit.theta <-x%*%beta1 
                                    logistic(logit.theta)
                                    } 
                # la proba de la reponse 
                dat2<-data.frame(dat,gm,gc);  
                outc1 <- rbinom(N, 1, apply(vx,1,thet1))
                datf <- data.frame(obs,outc1,dat2)
                nva1<-c("obs",nva)
                names(datf)<-nva1
                return(datf)
                }
                }
