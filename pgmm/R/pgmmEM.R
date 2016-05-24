.packageName<-'pgmm'
system.time({
m_best <-NA
G_best<-0
init_load<-function(x4,z4,G4,p4,q4){
    sampcov<-list(); 
    for(g in 1:G4){sampcov[[g]]<-matrix(0,p4,p4);}
    mu<-matrix(0,G4,p4)
    n<-colSums(z4)
    pi<-colSums(z4)/dim(x4)[1]
    for(g in 1:G4){mu[g,]<-apply(x4*z4[,g],2,sum)/(n[g]);}
    for(g in 1:G4){sampcov[[g]]<-cov.wt(x4,wt=(z4[,g]),center=mu[g,])$cov*(sum(z4[,g]^2)-1)/sum(z4[,g]);}
    lambda<-list();
    lambda1_tmp<-c(rep(0,p4*q4*G4));
    s<-1;
    for(g in 1:G4){
        evec<-eigen(sampcov[[g]])$vector
        eval<-eigen(sampcov[[g]])$values
        for(i in 1:p4){
            for(j in 1:q4){
                lambda1_tmp[s]<-sqrt(eval[j])*evec[i];
                s<-s+1;
            }
        }
    }
    lambda[["sep"]]<-lambda1_tmp
    
    psi<-list()
    lam_mat<-list()
    k4<-0
    for(g in 1:G4){
        lam_mat[[g]]<-matrix(lambda1_tmp[(1+k4):(g*p4*q4)],nrow=p4,ncol=q4,byrow=TRUE)
        k4<-p4*q4*g
    }
    temp_p6<-c(rep(0,p4))
    for(g in 1:G4){
        temp_p6 <- temp_p6 + pi[g]*abs(diag(sampcov[[g]]-lam_mat[[g]]%*%t(lam_mat[[g]])))
    }
    psi[[6]]<-temp_p6
    psi[[5]]<-sum(psi[[6]])/p4
    psi_tmp<-matrix(0,G4,p4)
    for (g in 1:G4){
        psi_tmp[g,]<-abs(diag(sampcov[[g]]-lam_mat[[g]]%*%t(lam_mat[[g]])))
    } 
    psi[[7]]<-rowMeans(psi_tmp)
    psi[[8]]<-as.vector(t(psi_tmp))
    
    stilde<-matrix(0,p4,p4)
    for(g in 1:G4){
        stilde<-stilde+pi[g]*sampcov[[g]];
    }
    evec<-eigen(stilde)$vector
    eval<-eigen(stilde)$values
    lambda_tilde<-matrix(0,p4,q4)
    s<-1
    for(i in 1:p4){
        for(j in 1:q4){
            lambda1_tmp[s]<-sqrt(eval[j])*evec[i];
            s<-s+1
        }
    }
    lambda[["tilde"]]<-lambda1_tmp
    
    lam_mat[[1]]<-matrix(lambda1_tmp[1:(p4*q4)],nrow=p4,ncol=q4,byrow=TRUE)
    psi[[2]]<-abs(diag(stilde-lam_mat[[1]]%*%t(lam_mat[[1]])))
    psi[[1]]<-sum(psi[[2]])/p4
    psi_tmp<-matrix(0,G4,p4)
    for (g in 1:G4){
        psi_tmp[g,]<-abs(diag(sampcov[[g]]-lam_mat[[1]]%*%t(lam_mat[[1]])))
    } 
    psi[[3]]<-rowMeans(psi_tmp)
    psi[[4]]<-as.vector(t(psi_tmp))
    psi[[9]]<-c(psi[[3]],rep(0,p4))
    psi[[10]]<-c(psi[[7]],rep(0,p4))
    psi[[11]]<-c(psi[[1]],rep(0,(G4*p4)))
    psi[[12]]<-c(psi[[5]],rep(0,(G4*p4)))
    
    lambda[["psi"]]<-psi
    lambda
}
endPrint<-function(icl,zstart,loop,m_best,q_best,G_best,bic_best,class_ind){
    start_names<-c("NA","k-means","custom")
    if(!class_ind){
        if(!icl){
            if(zstart==1){
                if(loop==1){
                    cat("Based on 1 random start, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, 
                        " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
                }else{
                    cat("Based on ", loop, " random starts, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
                }
            }else{
                cat("Based on ", start_names[zstart], " starting values, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
            }
        }else{
            if(zstart==1){
                if(loop==1){
                    cat("Based on 1 random start, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, 
                        " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
                }else{
                    cat("Based on ", loop, " random starts, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
                }
            }else{
                cat("Based on ", start_names[zstart], " starting values, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
            }
        }
    }else{
        if(!icl){
            cat("Based on the labelled and unlabelled data provided, the best model (BIC) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe BIC for this model is ", bic_best, ".", sep="")
        }else{
            cat("Based on the labelled and unlabelled data provided, the best model (ICL) for the range of factors and components used is a ", m_best, " model with q = ", q_best, " and G = ", G_best, ".\nThe ICL for this model is ", bic_best, ".", sep="")
        }
    }        
}
pgmmEM<-function(x,rG=1:2,rq=1:2,class=NULL,icl=FALSE,zstart=2,cccStart=TRUE,loop=3,zlist=NULL,modelSubset=NULL,seed=123456,tol=0.1,relax=FALSE){
	set.seed(seed)
	x<-as.matrix(x)
	run_pgmm<-function(x,z,bic,cls,q8,p,G8,N,model,cluster,lambda,psi,TOL=tol){
		p4<-.C("pgmm_c",as.double(x),as.double(z),as.double(bic),as.integer(cls),as.integer(q8),as.integer(p),as.integer(G8),as.integer(N),
			   as.integer(model),as.integer(cluster),as.double(lambda),as.double(psi),as.double(TOL),PACKAGE="pgmm")
		list(p4[[2]],p4[[3]],p4[[11]],p4[[12]])
	}
	is.int<-function(no){
		abs(no-round(no))<1e-15
	}
	if((tol>0.1)||(!is.double(tol))){stop("Invalid entry for tol.")}
	models_all<-c("CCC","CCU","CUC","CUU","UCC","UCU","UUC","UUU","CCUU","UCUU","CUCU","UUCU")
	model_num<-1:12
	names(model_num)<-models_all
	if(is.null(modelSubset)){modelSubset<-models_all}
	bic_out<-list()
	if(!is.null(zlist)||zstart==3){if(!is.list(zlist)){stop("Expected a list for zlist.")}}
	if(!is.null(class)){if(any(!is.double(class))&any(!is.int(class))){stop("The vector class may contain integers only.")}}
	Gmin<-rG[1]
	Gmax<-rG[length(rG)]
	qmin<-rq[1]
	qmax<-rq[length(rq)]
	G_offset<-Gmin-1
	q_offset<-qmin-1
	N<-dim(x)[1]
	p<- dim(x)[2]
	x1<-as.vector(t(x))
	bic_max<--Inf
	bic_best<--Inf
	if(!is.null(zlist)){
		for(g1 in rG){
			if(g1>1){
				if((any(!is.integer(zlist[[g1]])))){stop("Each element of zlist (G>1) must contain integers only.");}
				if(length(zlist[[g1]])!=N){stop("Each element of zlist (G>1) must have length equal to the number of samples.");}
			}
		}
	}
	test_paras<-c(qmax,qmin,Gmax,Gmin)
	if(!is.null(class))if(length(class)!=N){stop("The vector class must have length equal to the number of samples.");}
	if(any(!is.double(test_paras))&any(!is.integer(test_paras))){stop("The parameters rG and rq take ranges of integer values only.")}
	if(Gmax<Gmin){stop("The first element in the range rG is larger than the last.");}
	if(qmax<qmin){stop("The first element in the range rq is larger than the last.");}
	if(!relax){
        if(qmax>p/2){stop("qmax>p/2 is not allowed; set relax=TRUE to relax this constraint.");}
    }else{
        if(qmax>(3*p/4)){stop("qmax>3p/4 is not allowed.");}
    }
	if(!is.double(x)){stop("All elements of x must have type double.");}
	if(!is.logical(cccStart)){stop("cccStart takes logical values only.");}
	if(is.null(class)){class<-rep(0,N);class_ind<-0;}else{class_ind<-1;}
	for(mod in modelSubset){
		bic_temp<-matrix(NA,Gmax-Gmin+1,qmax-qmin+1)
		rownames(bic_temp)<-c(rG)
		colnames(bic_temp)<-c(rq)	
		bic_out[[mod]]<-bic_temp
	}
	if(class_ind){
        if(max(class)>Gmin){stop("The lowest value in rG cannot be less than max(class).")}
		for(g1 in rG){
			zt<-matrix(0,N,g1)
			cls_ind<-(class==0)
			for (i in 1:N){
				if(cls_ind[i]){zt[i,]<-1/g1}
				else{zt[i,class[i]]<-1}	
			}
			z1<-as.vector(t(zt))
			LAMBDA<-list()
			for (q1 in rq){
				LAMBDA[[q1-q_offset]]<-init_load(x,zt,g1,p,q1)
			}		
			for (m in modelSubset){
				for (q1 in rq){
                    if(substr(m,1,1)=="C"){
                        lam_temp<-LAMBDA[[q1-q_offset]]$tilde
                    }else{
                        lam_temp<-LAMBDA[[q1-q_offset]]$sep
                    }
                    psi_temp<-LAMBDA[[q1-q_offset]]$psi[[model_num[m]]]
                    temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,model_num[m],class_ind,lam_temp,psi_temp)
					bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
					if(!is.nan(temp[[2]])){
						if(icl&(g1>1)){
                            z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
                            mapZ<-rep(0,N)
                            for(i9 in 1:N){
                                mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
                            }
                            icl2<-0
                            for(i9 in 1:N){
                                icl2<-icl2+log(z_mat_tmp[i9,mapZ[i9]])
                            }
                            bic_out[[m]][g1-G_offset,q1-q_offset]<-bic_out[[m]][g1-G_offset,q1-q_offset]+icl2
                        }
                        if(temp[[2]]>bic_max){
							z_best<-temp[[1]];bic_best<-temp[[2]];
							bic_max<-bic_best;G_best<-g1;q_best<-q1;
							m_best<-m;lambda_best<-temp[[3]];psi_best<-temp[[4]]
						}
					}
				}
			}
		}
	}else{
		bic_start<-matrix(NA,Gmax-Gmin+1,qmax-qmin+1)
		if(zstart==1){
			for(l in 1:loop){
				for(g1 in rG){
					z<-matrix(0,N,g1)
                    for(i in 1:N){
                        sum<-0
                        for (j in 1:g1){
                            z[i,j]<-runif(1,0,1);
                            sum<-sum+z[i,j];
                        }
                        for (j in 1:g1){
                            z[i,j]<-z[i,j]/sum
                        }
                    }
					z1<-as.vector(t(z))
					LAMBDA<-list()
					for (q1 in rq){
						LAMBDA[[q1-q_offset]]<-init_load(x,z,g1,p,q1)
					}
					if(cccStart){
						bic_ccc_max<--Inf
						for (q1 in rq){
							temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,1,class_ind,LAMBDA[[q1-q_offset]]$tilde,LAMBDA[[q1-q_offset]]$psi[[1]])
							bic_start[g1-G_offset,q1-q_offset]<-temp[[2]]
							if(!is.nan(temp[[2]])){
								if(icl&(g1>1)){
									z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
									mapZ<-rep(0,N)
									for(i9 in 1:N){
										mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
									}
									icl1<-0
									if(icl){
                                        for(i9 in 1:N){
                                            icl1<-icl1+log(z_mat_tmp[i9,mapZ[i9]])
                                        }
                                    }
 									bic_start[g1-G_offset,q1-q_offset]<-bic_start[g1-G_offset,q1-q_offset]+icl1
								}
								if(bic_start[g1-G_offset,q1-q_offset]>bic_ccc_max){
									z_init_best<-temp[[1]];
                                    bic_ccc_max<-bic_start[g1-G_offset,q1-q_offset]
								}
							}
						}
						z_init_mat<-matrix(z_init_best,nrow=N,ncol=g1,byrow=TRUE);						
						for (q1 in rq){
							LAMBDA[[q1-q_offset]]<-init_load(x,z_init_mat,g1,p,q1)
						}
						z1<-as.vector(t(z_init_mat))
					}
                    for (m in modelSubset){
                        for (q1 in rq){
                            if(substr(m,1,1)=="C"){
                                lam_temp<-LAMBDA[[q1-q_offset]]$tilde
                            }else{
                                lam_temp<-LAMBDA[[q1-q_offset]]$sep
                            }
                            psi_temp<-LAMBDA[[q1-q_offset]]$psi[[model_num[m]]]
                            temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,model_num[m],class_ind,lam_temp,psi_temp)
                            if(!is.nan(temp[[2]])){
                            	if(icl&(g1>1)){
                                    z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
                                    mapZ<-rep(0,N)
                                    for(i9 in 1:N){
                                        mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
                                    }
                                    icl2<-0
                                    for(i9 in 1:N){
                                        icl2<-icl2+log(z_mat_tmp[i9,mapZ[i9]])
                                    }
                                    temp[[2]]<-temp[[2]]+icl2
                                }
                                if(l==1){
                            		bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
                            	}else if(is.na(bic_out[[m]][g1-G_offset,q1-q_offset])){
                            		bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
                            	}else if(bic_out[[m]][g1-G_offset,q1-q_offset]<temp[[2]]){
                            		bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
                            	}
                                if(temp[[2]]>bic_max){
                                    z_best<-temp[[1]];bic_best<-temp[[2]];
                                    bic_max<-bic_best;G_best<-g1;q_best<-q1;
                                    m_best<-m;lambda_best<-temp[[3]];psi_best<-temp[[4]]
                                }
                            }
                        }
                    }
				}
			}
		}else if((zstart==2)||(zstart==3)){
			for(g1 in rG){
				z<-matrix(0,N,g1)
				if(zstart==3){
					if(g1==1){z_ind<-c(rep(1,N))}
					else {z_ind<-zlist[[g1]]}
				}
				if(zstart==2){
					if(g1==1){z_ind<-c(rep(1,N))}
					else{set.seed(123456);z_ind<-kmeans(x,g1,nstart=5)$cluster;}
				}
				for (i in 1:N){
					z[i,z_ind[i]]<-1	
				}
				LAMBDA<-list()
				for (q1 in rq){
					LAMBDA[[q1-q_offset]]<-init_load(x,z,g1,p,q1)
				}
				z1<-as.vector(t(z))
				for (m in modelSubset){
					for (q1 in rq){
						if(substr(m,1,1)=="C"){
							lam_temp<-LAMBDA[[q1-q_offset]]$tilde
						}else{
							lam_temp<-LAMBDA[[q1-q_offset]]$sep
						}
						psi_temp<-LAMBDA[[q1-q_offset]]$psi[[model_num[m]]]
                        temp<-run_pgmm(x1,z1,0,class,q1,p,g1,N,model_num[m],class_ind,lam_temp,psi_temp)
						bic_out[[m]][g1-G_offset,q1-q_offset]<-temp[[2]]
						if(!is.nan(temp[[2]])){
							if(icl&(g1>1)){
								z_mat_tmp<-matrix(temp[[1]],nrow=N,ncol=g1,byrow=TRUE)
								mapZ<-rep(0,N)
								for(i9 in 1:N){
									mapZ[i9]<-which(z_mat_tmp[i9,1:g1]==max(z_mat_tmp[i9,1:g1]))
								}
								icl2<-0
								for(i9 in 1:N){
									icl2<-icl2+log(z_mat_tmp[i9,mapZ[i9]])
								}
								bic_out[[m]][g1-G_offset,q1-q_offset]<-bic_out[[m]][g1-G_offset,q1-q_offset]+icl2
							}
							if(temp[[2]]>bic_max){
								z_best<-temp[[1]];bic_best<-bic_out[[m]][g1-G_offset,q1-q_offset];
#	                                                        z_mat<-matrix(z_best,nrow=N,ncol=G_best,byrow=TRUE)
								bic_max<-bic_best;G_best<-g1;q_best<-q1;
	                                                        z_mat<-matrix(z_best,nrow=N,ncol=G_best,byrow=TRUE)
								m_best<-m;lambda_best<-temp[[3]];psi_best<-temp[[4]]
							}
						}else{z_mat<-NULL}
                                             
					}
				}
			}
		}else{
            stop("Invalid entry for zstart: 1 random; 2 k-means; 3 user-specified list.")
        }
	}
#	z_mat<-matrix(z_best,nrow=N,ncol=G_best,byrow=TRUE)
#        if((!is.nan(temp[[2]]))&temp[[2]]>bic_max){ 
        if((!is.nan(temp[[2]]))){ 
	if(substr(m_best,1,1)=="C"){
        lambda_mat<-matrix(lambda_best,nrow=p,ncol=q_best,byrow=TRUE)
    }else{
        lambda_mat<-NULL;
        for(g1 in 1:G_best){
            upper<-(q_best*p)*g1;
            lambda_mat[[g1]]<-matrix(lambda_best[(upper-(p*q_best-1)):upper],nrow=p,ncol=q_best,byrow=TRUE);
        }
    }
    if((m_best=="CUU")||(m_best=="UUU")){
        psi_mat<-list()
        for(g1 in 1:G_best){
            upper<-p*g1;
            psi_mat[[g1]]<-diag(psi_best[(upper-p+1):upper]);
        }
    }else if((m_best=="CCC")||(m_best=="UCC")){
        psi_mat<-psi_best[1];
    }else if((m_best=="CCU")||(m_best=="UCU")){
        psi_mat<-psi_best[1:p];
    }else if((m_best=="CUC")||(m_best=="UUC")){
        psi_mat<-list()
        for(g1 in 1:G_best){
            psi_mat[[g1]]<-psi_best[g1];
        }
    }else if((m_best=="CCUU")||(m_best=="UCUU")){
        psi_mat<-list();
        psi_mat[["omega"]]<-psi_best[1:G_best];
        psi_mat[["delta"]]<-diag(psi_best[(G_best+1):(G_best+p)]);
    }else if((m_best=="CUCU")||(m_best=="UUCU")){
        psi_mat<-list();
        psi_mat[["omega"]]<-psi_best[1];
        for(g1 in 1:G_best){
            temp_string<-paste("delta",toString(g1),sep="");
            lower<-2+(g1-1)*p
            psi_mat[[temp_string]]<-diag(psi_best[lower:(lower+p-1)]);
        }
    }
    z_mat<-matrix(z_best,nrow=N,ncol=G_best,byrow=TRUE)
  }else{z_mat=NULL  
}
    if(G_best >0){
    class_best<-rep(0,N)
	for(i in 1:N){
		class_best[i]<-which(z_mat[i,1:G_best]==max(z_mat[i,1:G_best]))
	}
    endPrint(icl,zstart,loop,m_best,q_best,G_best,bic_best,class_ind)
	if(!icl){
		foo<-list(map=class_best,model=m_best,g=G_best,q=q_best,bic=bic_out,zhat=z_mat,load=lambda_mat,noisev=psi_mat,plot_info=list(Gmin,Gmax,modelSubset,icl),summ_info=list(icl,zstart,loop,bic_best,class_ind))
	}else{
		foo<-list(map=class_best,model=m_best,g=G_best,q=q_best,icl=bic_out,zhat=z_mat,load=lambda_mat,noisev=psi_mat,plot_info=list(Gmin,Gmax,modelSubset,icl),summ_info=list(icl,zstart,loop,bic_best,class_ind))
	}
    class(foo)<-"pgmm"
    foo
} else{stop}
  }
summary.pgmm<-function(object,...){
    a<-object$summ_info[[1]]
    b<-object$summ_info[[2]]
    c<-object$summ_info[[3]]
    d<-object$model
    e<-object$q
    f<-object$g
    k<-object$summ_info[[4]]
    m<-object$summ_info[[5]]
    endPrint(a,b,c,d,e,f,k,m)
}
print.pgmm<-function(x,...){
    bicl<-x$plot_info[[4]]
    if(!bicl){
        cat("BIC for each model, number of components (rows), and number of factors (columns).\n")
        print.default(x$bic)
    }else{
        cat("ICL for each model, number of components (rows), and number of factors (columns).\n")
        print.default(x$icl)
    }
}
plot.pgmm<-function(x,onlyAll=FALSE,...){
    x$plot_info[[3]]->models
    x$plot_info[[4]]->icl1
	if(length(models)<3){
		par(mfrow=c(1,2),ask=FALSE)
	}else if(length(models)<5){
		par(mfrow=c(2,2),ask=FALSE)
	}else if(length(models)<7){
		par(mfrow=c(2,3),ask=FALSE)
	}else if(length(models)<10){
		par(mfrow=c(3,3),ask=FALSE)
	}else{
		par(mfrow=c(3,4),ask=FALSE)
	}
	if(icl1){
        bicl<-x$icl
        ylabel<-"ICL"
    }else{
        bicl<-x$bic
        ylabel<-"BIC"
    }
    for(k in models){
        n<-which(models==k)
        matplot(c(x$plot_info[[1]]:x$plot_info[[2]]),bicl[[n]],type="b",ylab=ylabel,xlab="G",main=k,xaxt="n")
        axis(1,at=c(x$plot_info[[1]]:x$plot_info[[2]]))
	}
    
    if(!onlyAll){
		for(k in models){
	        par(mfrow=c(1,1),ask=TRUE)
	        n<-which(models==k)
    	    matplot(c(x$plot_info[[1]]:x$plot_info[[2]]),bicl[[n]],type="b",ylab=ylabel,xlab="G",main=k,xaxt="n")
        	axis(1,at=c(x$plot_info[[1]]:x$plot_info[[2]]))
        }
	}
}
})
