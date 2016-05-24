# Version 9/22/2015
upsilon <-function(data,x,m,y,covs=NULL,parallel.med=NULL,seq.med=NULL)
{
if(!requireNamespace("lavaan", quietly = TRUE)) stop("The package 'lavaan' is needed; please install the package and try again.")
  
  UpsES <- list()
     class(UpsES) <- 'UpsES'
     if(length(m)==1){ #single mediator model
          
          ####Total indirect effects and upsilon####
          UpsES$IE$Total <- list()
          UpsES$Upsilon$Total <- list()
          bx <- as.numeric(length(x))
          
          ###create lavaan model input
          modely <- paste0(c(y,paste0(c(m,paste0(x,collapse='+'),covs),collapse='+')),collapse='~')
          modelm <- paste0(c(m,paste0(c(paste0(x,collapse='+'),covs),collapse='+')),collapse='~')
          model <- paste0(c(modely,modelm),collapse='\n ')
          
          ###compute coefs in lavaan
          out <- lavaan::sem(model,data=data)
          
          ###covariance matrix from lavaan output
          covmat <- lavaan::inspect(out,'cor.all')
          
          ###B matrix 
          B <- lavaan::inspect(out,'coef')$beta
          B <- unclass(B)
          I <- diag(dim(B)[1])
          
          ###Indirect effect matrix 
          Ind <- solve(I-B)-I-B
          
          ###Indirect effects
          for(i in seq(bx)){
               UpsES$IE$Total[[i]] <- Ind[y,x[i]]
               names(UpsES$IE$Total[[i]]) <- paste0(c(x[i],m,y),collapse="->")
               if(!is.null(covs)){
                    names(UpsES$IE$Total[[i]]) <- paste0(c(names(UpsES$IE$Total[i]),covs),collapse='|')
               }
          }
          
          ###upsilon matrix
          ups <- Ind%*%covmat%*%t(Ind)
          UpsES$Upsilon$Total[[1]] <- ups[y,y]
          names(UpsES$Upsilon$Total[[1]]) <- paste0(c(paste0(x,collapse='+'),m,y),collapse="->")
          if(!is.null(covs)){
               names(UpsES$Upsilon$Total[[1]]) <- paste0(c(names(UpsES$Upsilon$Total[[1]]),covs),collapse='|')
          }
          
          ####Specific and Unconditional SSIEs####
          if(bx>1){
               UpsES$IE$Specific <- list()
               UpsES$IE$Uncond <- list()
               UpsES$Upsilon$Specific <- list()
               UpsES$Upsilon$Uncond <- list()
               UpsES$Upsilon$Unique <- list()
               for(i in seq(bx)){
                    O <- diag(bx+length(covs)+2)
                    colnames(O) <- rownames(O) <- rownames(B)
                    O[,setdiff(c(x,m,y,covs),c(x[i],m,y))] <- 0
                    Bstar_sp <- t(O)%*%B%*%O
                    Indstar_sp <- solve(I-Bstar_sp)-I-Bstar_sp
                    
                    ####Specific indirect effets
                    UpsES$IE$Specific[[i]] <- Indstar_sp[y,x[i]]
                    names(UpsES$IE$Specific[[i]]) <- paste0(c(x[i],m,y),collapse="->")
                    if(!is.null(covs)){
                         names(UpsES$IE$Specific[[i]]) <- paste0(c(names(UpsES$IE$Specific[[i]]),covs),collapse='|')
                    }
                    
                    ###Specific SSIEs
                    ups_sp <- Indstar_sp%*%covmat%*%t(Indstar_sp)
                    UpsES$Upsilon$Specific[[i]] <- ups_sp[y,y]
                    names(UpsES$Upsilon$Specific[[i]]) <- paste0(c(x[i],m,y),collapse="->")
                    if(!is.null(covs)){
                         names(UpsES$Upsilon$Specific[[i]]) <- paste0(c(names(UpsES$Upsilon$Specific[[i]]),covs),collapse='|')
                    }
                    
                    ###unconditional models
                    modely_unc <- paste0(c(y,paste0(c(m,paste0(x[-i],collapse='+'),covs),collapse='+')),collapse='~')
                    modelm_unc <- paste0(c(m,paste0(c(paste0(x[-i],collapse='+'),covs),collapse='+')),collapse='~')
                    model_unc <- paste0(c(modely_unc,modelm_unc),collapse='\n ')
                    out_unc <- lavaan::sem(model_unc,data=data)
                    covmat_unc <- lavaan::inspect(out_unc,'cor.all')
                    B_unc <- lavaan::inspect(out_unc,'coef')$beta
                    B_unc <- unclass(B_unc)
                    I_unc <- diag(dim(B_unc)[1])
                    Ind_unc <- solve(I_unc-B_unc)-I_unc-B_unc
                    
                    #Unconditional indirect effets
                    ux <- x[-i]
                    for(j in seq(length(ux))){
                         UpsES$IE$Uncond[[j+(i-1)*length(ux)]] <- Ind_unc[y,ux[j]]
                         names(UpsES$IE$Uncond[[j+(i-1)*length(ux)]]) <- paste0(c(paste0(c(paste0(ux[j],collapse='+'),m,y),collapse="->"),x[i]),collapse=' w/o ')
                         if(!is.null(covs)){
                              names(UpsES$IE$Uncond[[j+(i-1)*length(ux)]]) <- paste0(c(names(UpsES$IE$Uncond[[j+(i-1)*length(ux)]]),covs),collapse='|')
                         }
                    }
                    
                    #Unconditional SSIEs
                    ups_unc <- Ind_unc%*%covmat_unc%*%t(Ind_unc)
                    UpsES$Upsilon$Uncond[[i]] <- ups_unc[y,y]
                    names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(paste0(c(paste0(x[-i],collapse='+'),m,y),collapse="->"),x[i]),collapse=' w/o ')
                    if(!is.null(covs)){
                         names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(names(UpsES$Upsilon$Uncond[[i]]),covs),collapse='|')
                    }
               }
               
               ####Unique SSIEs####
               for(i in seq(bx)){
                    UpsES$Upsilon$Unique[[i]] <- UpsES$Upsilon$Total[[1]] - UpsES$Upsilon$Uncond[[i]]
                    parse <- unlist(strsplit(names(UpsES$Upsilon$Uncond[[i]]),'->'))
                    if(length(grep('\\+',parse,value=TRUE)!=0)){
                         parse <- grep('\\+',parse,value=TRUE)
                    }
                    names(UpsES$Upsilon$Unique[[i]]) <- paste0(c(x[which(!x %in% unlist(strsplit(parse,'+',fixed=TRUE)))],m,y),collapse="->")
                    if(!is.null(covs)){
                         names(UpsES$Upsilon$Unique[[i]]) <- paste0(c(names(UpsES$Upsilon$Unique[[i]]),covs),collapse='|')
                    }
               }
          }
     } else if(length(m)>1){ #i.e. mediation with multiple mediator
          bx <- as.numeric(length(x))
          bm <- as.numeric(length(m))
          
          UpsES$IE$Total <- list()
          UpsES$IE$Specific <- list()
          UpsES$IE$Uncond <- list()
          
          UpsES$Upsilon$Total <- list()
          UpsES$Upsilon$Specific <- list()
          UpsES$Upsilon$Uncond  <- list()
          UpsES$Upsilon$Unique <- list()
          
          
          ####Total indirect effects and upsilon
          ###create lavaan model input
          models <- list()
          usedmeds <- list()
          lav.models <- list()
          used <- NULL
          if(!is.null(seq.med)){
               for(i in seq(length(seq.med))){
                    seqmeds <- unlist(strsplit(unlist(strsplit(seq.med[[i]],'~',fixed=TRUE)),'+',fixed=TRUE))
                    if(!is.null(parallel.med)){     
                         if(any(seqmeds %in% unlist(strsplit(unlist(parallel.med),'+',fixed=TRUE)))){
                              for(j in seq(length(parallel.med))){
                                   parmeds <- unlist(strsplit(parallel.med[[j]],'+',fixed=TRUE))
                                   if(any(seqmeds %in% parmeds)){
                                        if(is.null(used)){
                                             if(length(which(seqmeds %in% parmeds))>1){
                                                  lav.models<- append(lav.models,paste0(c(seqmeds[!seqmeds %in% parmeds],paste0(c(seqmeds[seqmeds %in% parmeds],x,covs),collapse='+')),collapse='~'))
                                             } else {
                                                  lav.models <- append(lav.models,paste0(c(seqmeds[!seqmeds %in% parmeds],paste0(c(x,covs),collapse='+')),collapse='~'))
                                                  lav.models <- append(lav.models,paste0(c(seqmeds[seqmeds %in% parmeds],paste0(c(seqmeds[!seqmeds %in% parmeds],x,covs),collapse='+')),collapse='~'))
                                             }
                                        } else {
                                             if(length(which(seqmeds %in% parmeds))>1){
                                                  lav.models <- append(lav.models,paste0(c(seqmeds[!seqmeds %in% parmeds],paste0(c(paste0(seqmeds[seqmeds %in% parmeds],collapse='+'),used[!used %in% seqmeds & !used %in% parallel.med[[j]]],x,covs),collapse='+')),collapse='~'))
                                             } else if(seqmeds[!seqmeds %in% parmeds] %in% used){
                                                  lav.models <- append(lav.models,paste0(c(seqmeds[seqmeds %in% parmeds],paste0(c(seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds[!seqmeds %in% parmeds]],x,covs),collapse='+')),collapse='~'))
                                             } else lav.models <- append(lav.models,paste0(c(seqmeds[seqmeds %in% parmeds],paste0(c(seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds],x,covs),collapse='+')),collapse='~'))
                                        }
                                        for(k in seq(length(x))){
                                             if(is.null(used)){
                                                  if(length(which(seqmeds %in% parmeds))>1){
                                                       models<- append(models,paste0(c(seqmeds[!seqmeds %in% parmeds],paste0(seqmeds[seqmeds %in% parmeds],collapse='+'),x[[k]]),collapse='<-'))
                                                  } else models <- append(models,paste0(c(seqmeds[seqmeds %in% parmeds],seqmeds[!seqmeds %in% parmeds],x[[k]]),collapse='<-'))
                                             } else {
                                                  if(length(which(seqmeds %in% parmeds))>1){
                                                       models <- append(models,paste0(c(seqmeds[!seqmeds %in% parmeds],paste0(seqmeds[seqmeds %in% parmeds],collapse='+'),used[!used %in% seqmeds & !used %in% parallel.med[[j]]],x[[k]]),collapse='<-'))
                                                  } else if(seqmeds[!seqmeds %in% parmeds] %in% used){
                                                       models <- append(models,paste0(c(seqmeds[seqmeds %in% parmeds],seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds[!seqmeds %in% parmeds]],x[[k]]),collapse='<-'))
                                                  } else models <- append(models,paste0(c(seqmeds[seqmeds %in% parmeds],seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds],x[[k]]),collapse='<-'))
                                             }
                                        }
                                        if(any(seqmeds %in% parmeds)){
                                             if(length(which(parmeds %in% seqmeds))==length(parmeds)){
                                                  usedmeds[[i]] <- rev(seqmeds)
                                                  usedmeds[[i]][which(seqmeds %in% parmeds)] <- parallel.med[[j]]
                                             } else usedmeds[[i]] <- seqmeds[!seqmeds %in% parmeds]
                                        } else usedmeds[[i]] <- seqmeds
                                   } 
                              }
                              
                         } else {
                              if(is.null(used)){
                                   lav.models <- append(lav.models,paste0(c(seqmeds[2],paste0(c(x,covs),collapse='+')),collapse='~'))
                                   lav.models <- append(lav.models,paste0(c(seqmeds[1],paste0(c(seqmeds[2],x,covs),collapse='+')),collapse='~'))
                              } else {
                                   lav.models <- append(lav.models,paste0(c(seqmeds[2],paste0(c(used[!used %in% seqmeds],x,covs),collapse='+')),collapse='~'))
                                   lav.models <- append(lav.models,paste0(c(seqmeds[1],paste0(c(seqmeds[2],used[!used %in% seqmeds],x,covs),collapse='+')),collapse='~'))
                              }
                              for(k in seq(length(x))){
                                   if(is.null(used)){
                                        models <- append(models,paste0(c(seqmeds[1],seqmeds[2],x[[k]]),collapse='<-'))
                                   } else models <- append(models,paste0(c(seqmeds[1],seqmeds[2],used[!used %in% seqmeds],x[[k]]),collapse='<-'))
                              }
                              usedmeds[[i]] <- seqmeds
                         } 
                    } else {
                         if(is.null(used)){
                              lav.models <- append(lav.models,paste0(c(seqmeds[2],paste0(c(x,covs),collapse='+')),collapse='~'))
                              lav.models <- append(lav.models,paste0(c(seqmeds[1],paste0(c(seqmeds[2],x,covs),collapse='+')),collapse='~'))
                         } else {
                              lav.models <- append(lav.models,paste0(c(seqmeds[1],paste0(c(seqmeds[2],used[!used %in% seqmeds],x,covs),collapse='+')),collapse='~'))
                         }
                         for(k in seq(length(x))){
                              if(is.null(used)){
                                   models <- append(models,paste0(c(seqmeds[1],seqmeds[2],x[[k]]),collapse='<-'))
                              } else models <- append(models,paste0(c(seqmeds[1],seqmeds[2],used[!used %in% seqmeds],x[[k]]),collapse='<-'))
                         }
                         usedmeds[[i]] <- rev(seqmeds)
                    }
                    used <- rev(unique(unlist(usedmeds),fromLast=TRUE))
               }
               if(!is.null(parallel.med)){
                    if(length(which(seqmeds %in% parmeds))>1){
                         lav.models <- append(lav.models,paste0(c(y,paste0(c(paste0(parmeds,collapse='+'),seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds & !used %in% parallel.med[[length(parallel.med)]]],x,covs),collapse='+')),collapse='~'))
                    } else if(any(seqmeds %in% parmeds)){
                         lav.models <- append(lav.models,paste0(c(y,paste0(c(parallel.med[[length(parallel.med)]],seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds[!seqmeds %in% parmeds]],x,covs),collapse='+')),collapse='~'))
                    } else lav.models <- append(lav.models,paste0(c(y,paste0(c(parallel.med[[length(parallel.med)]],seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds],x,covs),collapse='+')),collapse='~'))
               } else lav.models <- append(lav.models,paste0(c(y, paste0(c(used,x,covs),collapse='+')),collapse='~'))
               for(i in seq(length(x))){
                    if(!is.null(parallel.med)){
                         if(length(which(seqmeds %in% parmeds))>1){
                              models <- append(models,paste0(c(y,seqmeds[!seqmeds %in% parmeds],paste0(parmeds,collapse='+'),used[!used %in% seqmeds & !used %in% parallel.med[[length(parallel.med)]]],x[[i]]),collapse='<-'))
                         }else if(any(seqmeds %in% parmeds)){
                              models <- append(models,paste0(c(y,parallel.med[[length(parallel.med)]],seqmeds[!seqmeds %in% parmeds],used[!used %in% seqmeds[!seqmeds %in% parmeds]],x[[i]]),collapse='<-'))
                         } else models <- append(models,paste0(c(y, seqmeds[!seqmeds %in% used],used,x[[i]]),collapse='<-'))
                    } else models <- append(models,paste0(c(y, used,x[[i]]),collapse='<-'))
               }
          } else {
               for(i in seq(bm)){
                    lav.models <- append(lav.models, paste0(c(m[[i]],paste0(c(x,covs),collapse='+')),collapse='~'))
               }
               lav.models <- append(lav.models,paste0(c(y,paste0(c(parallel.med[[1]],x,covs),collapse='+')),collapse='~'))
               for(i in seq(length(x))){
                    models <- append(models,paste0(c(y,parallel.med[[1]],x[[i]]),collapse='<-'))
               }
          }
          
          lav.model <- paste0(rev(lav.models),collapse='\n ')
          
          ###compute coefs in lavaan
          out <- lavaan::sem(lav.model,data=data)
          
          ###covariance matrix from lavaan output
          covmat <- lavaan::inspect(out,'cor.all')
          ###B matrix from lavaan output, needs to be reordered
          B <- lavaan::inspect(out,'coef')$beta
          B <- unclass(B)
          O <- diag(bx+bm+length(covs)+1)
          colnames(O) <- rownames(O) <- rownames(B)
          O[,setdiff(c(x,m,y,covs),c(x,m,y))] <- 0
          I <- diag(dim(B)[1])
          Bstar_tot <- t(O)%*%B%*%O
          ###Indirect effect matrix 
          Ind_tot <- solve(I-Bstar_tot)-I-Bstar_tot
          
          tot_models <- matrix(0,nrow=length(models),ncol=2)
          for(i in seq(dim(tot_models)[1])){
               fx <- unlist(strsplit(models[[i]],'<-',fixed=TRUE))
               tot_models[i,] <- c(fx[1],fx[length(fx)]) 
          }
          for(i in seq(dim(tot_models)[1])){
               UpsES$IE$Total[[i]] <- Ind_tot[tot_models[i,1],tot_models[i,2]]
               names(UpsES$IE$Total[[i]]) <- paste0(rev(unlist(strsplit(models[[i]],'<-',fixed=TRUE))),collapse='->')
               if(!is.null(covs)){
                    names(UpsES$IE$Total[[i]]) <- paste0(c(names(UpsES$IE$Total[i]),covs),collapse='|')
               }
          }
          
          ups <- Ind_tot%*%covmat%*%t(Ind_tot)
          upstot_outs <- unique(tot_models[,1])
          for(i in seq(length(upstot_outs))){
               UpsES$Upsilon$Total[[i]] <- ups[upstot_outs[i],upstot_outs[i]]
               x.switch <- unlist(strsplit(models[[i*bx]],'<-',fixed=TRUE))
               names(UpsES$Upsilon$Total[[i]]) <- paste0(c(paste0(x,collapse='+'),rev(x.switch[!x.switch %in% x])),collapse='->')
               if(!is.null(covs)){
                    names(UpsES$Upsilon$Total[[i]]) <- paste0(c(names(UpsES$Upsilon$Total[i]),covs),collapse='|')
               }
          }
          
          if(length(seq.med)>0 & length(parallel.med)==0){
               iecount <- 0
               upscount <- 0
               for(i in seq(bx)){
                    for(j in seq(length(models))){
                         xmodel.sp <- unlist(strsplit(models[[j]],'<-',fixed=TRUE))
                         if(length(which(!xmodel.sp %in% x[[i]]))>2){
                              iecount <- iecount+1
                              upscount <- upscount+1
                              O <- diag(bx+bm+length(covs)+1)
                              colnames(O) <- rownames(O) <- rownames(B)
                              O[,setdiff(c(x,m,y),xmodel.sp[!xmodel.sp %in% x[[i]]])] <- 0
                              Bstar_sp <- O%*%B%*%O
                              I <- diag(dim(Bstar_sp)[1])
                              Indstar_sp <- solve(I-Bstar_sp)-I-Bstar_sp
                              if(length(which(xmodel.sp %in% x & !xmodel.sp %in% x[[i]]))>1){
                                   for(k in seq(length(which(!x %in% x[[i]])))){
                                        UpsES$IE$Specific[[iecount]] <- Indstar_sp[xmodel.sp[[1]],x[!x %in% x[[i]]]]
                                        names(UpsES$IE$Specific[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% x[[i]]]),collapse='->')
                                        if(!is.null(covs)){
                                             names(UpsES$IE$Specific[[iecount]]) <- paste0(c(names(UpsES$IE$Specific[iecount]),covs),collapse='|')
                                        }
                                        iecount<-iecount+1
                                   }
                              } else {
                                   UpsES$IE$Specific[[iecount]] <- Indstar_sp[xmodel.sp[[1]],xmodel.sp[[length(xmodel.sp[!xmodel.sp %in% x[[i]]])]]]
                                   names(UpsES$IE$Specific[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% x[[i]]]),collapse='->')
                                   if(!is.null(covs)){
                                        names(UpsES$IE$Specific[[iecount]]) <- paste0(c(names(UpsES$IE$Specific[iecount]),covs),collapse='|')
                                   }
                              }
                              ups_sp <- Indstar_sp%*%covmat%*%t(Indstar_sp)
                              UpsES$Upsilon$Specific[[upscount]] <- ups_sp[xmodel.sp[[1]],xmodel.sp[[1]]]
                              ups_sp.names <- rev(xmodel.sp[!xmodel.sp %in% x[[i]]])
                              if(length(bx)>1){
                                   names(UpsES$Upsilon$Specific[[upscount]]) <- paste0(c(paste0(x[x %in% ups_sp.names],collapse='+'),ups_sp.names[!ups_sp.names %in% x]),collapse='->')
                              } else names(UpsES$Upsilon$Specific[[upscount]]) <- paste0(ups_sp.names[!ups_sp.names %in% x],collapse='->')
                              if(!is.null(covs)){
                                   names(UpsES$Upsilon$Specific[[upscount]]) <- paste0(c(names(UpsES$Upsilon$Specific[upscount]),covs),collapse='|')
                              }
                              
                              umodels <- list()
                              for(k in seq(length(models))){
                                   umod <- unlist(strsplit(unlist(strsplit(models[[k]],'<-',fixed=TRUE)),'+',fixed=TRUE))
                                   umod <- umod[!umod %in% x[[i]]]
                                   if(length(which(umod %in% x))>0){
                                        xs <- x[!x %in% x[[i]]]
                                        for(l in seq(length(xs))){
                                             umodels <- append(umodels,paste0(c(umod,xs[[l]]),collapse='<-'))
                                        }
                                   } else {
                                        umodels <- append(umodels,paste0(umod,collapse='<-'))
                                   }
                              }                               
                              ulav.models <- list()
                              for(k in seq(length(lav.models))){
                                   ulmod <- unlist(strsplit(unlist(strsplit(lav.models[[k]],'~',fixed=TRUE)),'+',fixed=TRUE))
                                   ulmod <- ulmod[!ulmod %in% x[[i]]]
                                   if(length(ulmod)>1){
                                        if(!any(paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~')==ulav.models)){
                                             ulav.models <- append(ulav.models,paste0(c(ulmod[1],paste0(c(ulmod[ulmod!=ulmod[1]],covs),collapse='+')),collapse='~'))
                                        }
                                   }
                              }
                              ulav.model <- paste0(rev(ulav.models),collapse='\n ')
                              out <- lavaan::sem(ulav.model,data=data)
                              ###covariance matrix from lavaan output
                              covmatunc <- lavaan::inspect(out,'cor.all')
                              Bunc <- lavaan::inspect(out,'coef')$beta
                              Bunc <- unclass(Bunc)
                              O <- diag(dim(Bunc)[1])
                              colnames(O) <- rownames(O) <- rownames(Bunc)
                              O[,setdiff(c(x,m,y,covs),c(x,m,y))] <- 0
                              Iunc <- diag(dim(Bunc)[1])
                              Bstar_unc <- t(O)%*%Bunc%*%O
                              
                              ###Indirect effect matrix 
                              Indstar_unc <- solve(Iunc-Bstar_unc)-Iunc-Bstar_unc
                              
                              if(length(which(xmodel.sp %in% x & !xmodel.sp %in% x[[i]]))>1){
                                   for(k in seq(length(which(!x %in% x[[i]])))){
                                        UpsES$IE$Uncond[[iecount]] <- Indstar_unc[xmodel.sp[[1]],x[!x %in% x[[i]]]]
                                        names(UpsES$IE$Uncond[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% x[[i]]]),collapse='->')
                                        if(!is.null(covs)){
                                             names(UpsES$IE$Uncond[[iecount]]) <- paste0(c(names(UpsES$Uncond$Specific[iecount]),covs),collapse='|')
                                        }
                                        iecount<-iecount+1
                                   }
                              } else {
                                   UpsES$IE$Uncond[[iecount]] <- Indstar_unc[xmodel.sp[[1]],xmodel.sp[[length(xmodel.sp[!xmodel.sp %in% x[[i]]])]]]
                                   names(UpsES$IE$Uncond[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% x[[i]]]),collapse='->')
                                   if(!is.null(covs)){
                                        names(UpsES$IE$Uncond[[iecount]]) <- paste0(c(names(UpsES$IE$Uncond[iecount]),covs),collapse='|')
                                   }
                              }
                              ups_unc <- Indstar_unc%*%covmatunc%*%t(Indstar_unc)
                              UpsES$Upsilon$Uncond[[upscount]] <- ups_unc[xmodel.sp[[1]],xmodel.sp[[1]]]
                              ups_unc.names <- rev(xmodel.sp[!xmodel.sp %in% x[[i]]])
                              if(length(bx)>1){
                                   names(UpsES$Upsilon$Uncond[[upscount]]) <- paste0(c(paste0(x[x %in% ups_unc.names],collapse='+'),ups_unc.names[!ups_unc.names %in% x]),collapse='->')
                              } else names(UpsES$Upsilon$Uncond[[upscount]]) <- paste0(ups_unc.names[!ups_unc.names %in% x],collapse='->')
                              if(!is.null(covs)){
                                   names(UpsES$Upsilon$Uncond[[upscount]]) <- paste0(c(names(UpsES$Upsilon$Uncond[upscount+(i-1)*upscount]),covs),collapse='|')
                              }
                         }
                    }
                    liesp <- length(UpsES$IE$Specific)
                    lusp <- length(UpsES$Upsilon$Specific)
                    for(i in seq(bm)){
                         for(j in seq(length(models))){
                              xmodel.sp <- unlist(strsplit(models[[j]],'<-',fixed=TRUE))
                              if(length(which(!xmodel.sp %in% m[[i]]))>2){
                                   if(length(xmodel.sp[!xmodel.sp %in% m[[i]]])!=length(xmodel.sp)){
                                        iecount <- iecount+1
                                        upscount <- upscount+1
                                        O <- diag(bx+bm+length(covs)+1)
                                        colnames(O) <- rownames(O) <- rownames(B)
                                        O[,setdiff(c(x,m,y),xmodel.sp[!xmodel.sp %in% m[[i]]])] <- 0
                                        Bstar_sp <- O%*%B%*%O
                                        I <- diag(dim(Bstar_sp)[1])
                                        Indstar_sp <- solve(I-Bstar_sp)-I-Bstar_sp
                                        
                                        for(k in seq(bx)){
                                             UpsES$IE$Specific[[iecount]] <- Indstar_sp[xmodel.sp[!xmodel.sp %in% m[[i]]][1],x[[k]]]
                                             names(UpsES$IE$Specific[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% m[[i]]]),collapse='->')
                                        }
                                        if(!is.null(covs)){
                                             names(UpsES$IE$Specific[[iecount]]) <- paste0(c(names(UpsES$IE$Specific[iecount]),covs),collapse='|')
                                        }
                                        
                                        ups_sp <- Indstar_sp%*%covmat%*%t(Indstar_sp)
                                        UpsES$Upsilon$Specific[[upscount]] <- ups_sp[xmodel.sp[!xmodel.sp %in% m[[i]]][1],xmodel.sp[!xmodel.sp %in% m[[i]]][1]]
                                        ups_sp.names <- rev(xmodel.sp[!xmodel.sp %in% m[[i]]])
                                        names(UpsES$Upsilon$Specific[[upscount]]) <- paste0(c(paste0(x,collapse='+'),ups_sp.names[ups_sp.names %in% m],ups_sp.names[!ups_sp.names %in% m & !ups_sp.names %in% x]),collapse='->')
                                        if(!is.null(covs)){
                                             names(UpsES$Upsilon$Specific[[upscount]]) <- paste0(c(names(UpsES$Upsilon$Specific[upscount]),covs),collapse='|')
                                        }
                                        umodels <- list()
                                        for(k in seq(length(models))){
                                             umod <- unlist(strsplit(unlist(strsplit(models[[k]],'<-',fixed=TRUE)),'+',fixed=TRUE))
                                             umod <- umod[!umod %in% m[[i]]]
                                             umodels <- append(umodels,paste0(c(umod),collapse='<-'))
                                        }
                                        ulav.models <- list()
                                        for(k in seq(length(lav.models))){
                                             ulmod <- unlist(strsplit(unlist(strsplit(lav.models[[k]],'~',fixed=TRUE)),'+',fixed=TRUE))
                                             ulmod <- ulmod[!ulmod %in% m[[i]]]
                                             if(length(ulmod)>1){
                                                  if(!any(paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~')==ulav.models)){
                                                       ulav.models <- append(ulav.models,paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~'))
                                                  }
                                             }
                                        }
                                        ulav.model <- paste0(rev(ulav.models),collapse='\n ')
                                        out <- lavaan::sem(ulav.model,data=data)
                                        ###covariance matrix from lavaan output
                                        covmatunc <- lavaan::inspect(out,'cor.all')
                                        Bunc <- lavaan::inspect(out,'coef')$beta
                                        Bunc <- unclass(Bunc)
                                        O <- diag(dim(Bunc)[1])
                                        colnames(O) <- rownames(O) <- rownames(Bunc)
                                        O[,setdiff(c(x,m,y,covs),c(x,m,y))] <- 0
                                        Iunc <- diag(dim(Bunc)[1])
                                        Bstar_unc <- t(O)%*%Bunc%*%O
                                        
                                        ###Indirect effect matrix 
                                        Indstar_unc <- solve(Iunc-Bstar_unc)-Iunc-Bstar_unc
                                        if(bx>1){
                                             for(k in seq(bx)){
                                                  UpsES$IE$Uncond[[iecount]] <- Indstar_unc[xmodel.sp[!xmodel.sp %in% m[[i]]][1],x[k]]
                                                  names(UpsES$IE$Uncond[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% m[[i]]]),collapse='->')
                                                  if(!is.null(covs)){
                                                       names(UpsES$IE$Uncond[[iecount]]) <- paste0(c(names(UpsES$Uncond$Specific[iecount]),covs),collapse='|')
                                                  }
                                                  iecount <- iecount+1
                                             } 
                                        } else {
                                             UpsES$IE$Uncond[[iecount]] <- Indstar_unc[xmodel.sp[!xmodel.sp %in% m[[i]]][1],x]
                                             names(UpsES$IE$Uncond[[iecount]]) <- paste0(rev(xmodel.sp[!xmodel.sp %in% m[[i]]]),collapse='->')
                                        }
                                        ups_unc <- Indstar_unc%*%covmatunc%*%t(Indstar_unc)
                                        UpsES$Upsilon$Uncond[[upscount]] <- ups_unc[xmodel.sp[!xmodel.sp %in% m[[i]]][1],xmodel.sp[!xmodel.sp %in% m[[i]]][1]]
                                        ups_unc.names <- rev(xmodel.sp[!xmodel.sp %in% m[[i]]])
                                        if(length(bx)>1){
                                             names(UpsES$Upsilon$Uncond[[upscount]]) <- paste0(c(paste0(x,collapse='+'),ups_unc.names[!ups_unc.names %in% x]),collapse='->')
                                        } else names(UpsES$Upsilon$Uncond[[upscount]]) <- paste0(ups_unc.names,collapse='->')
                                        if(!is.null(covs)){
                                             names(UpsES$Upsilon$Uncond[[upscount]]) <- paste0(c(names(UpsES$Upsilon$Uncond[upscount]),covs),collapse='|')
                                        }
                                   }
                              }
                         }
                    }
               }
          } else if(length(seq.med)==0 & length(parallel.med)>0){ #parallel mediators
               if(bx>1){
                    for(i in seq(bx)){
                         xmodel.sp <- unlist(strsplit(unlist(strsplit(models[[i]],'<-',fixed=TRUE)),'+',fixed=TRUE))
                         O <- diag(bx+bm+length(covs)+1)
                         colnames(O) <- rownames(O) <- rownames(B)
                         O[,setdiff(c(x,m,y),c(x[[i]],m,y))] <- 0
                         Bstar_sp <- O%*%B%*%O
                         I <- diag(dim(Bstar_sp)[1])
                         Indstar_sp <- solve(I-Bstar_sp)-I-Bstar_sp
                         UpsES$IE$Specific[[i]] <- Indstar_sp[xmodel.sp[[1]],x[i]]
                         names(UpsES$IE$Specific[[i]]) <- paste0(c(x[[i]],parallel.med[[1]],xmodel.sp[[1]]),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$IE$Specific[[i]]) <- paste0(c(names(UpsES$IE$Specific[i]),covs),collapse='|')
                         }
                         ups_sp <- Indstar_sp%*%covmat%*%t(Indstar_sp)
                         UpsES$Upsilon$Specific[[i]] <- ups_sp[xmodel.sp[[1]],xmodel.sp[[1]]]
                         names(UpsES$Upsilon$Specific[[i]]) <- paste0(c(x[[i]],parallel.med[[1]],y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$Upsilon$Specific[[upscount]]) <- paste0(c(names(UpsES$Upsilon$Specific[upscount]),covs),collapse='|')
                         }
                         
                         umodels <- models[!models %in% models[i]]
                         
                         ulav.models <- list()
                         for(k in seq(length(lav.models))){
                              ulmod <- unlist(strsplit(unlist(strsplit(lav.models[[k]],'~',fixed=TRUE)),'+',fixed=TRUE))
                              ulmod <- ulmod[!ulmod %in% x[[i]]]
                              ulav.models <- append(ulav.models,paste0(c(ulmod[1],paste0(c(ulmod[ulmod!=ulmod[1]],covs),collapse='+')),collapse='~'))
                         }
                         ulav.model <- paste0(rev(ulav.models),collapse='\n ')
                         out <- lavaan::sem(ulav.model,data=data)
                         ###covariance matrix from lavaan output
                         covmatunc <- lavaan::inspect(out,'cor.all')
                         Bunc <- lavaan::inspect(out,'coef')$beta
                         Bunc <- unclass(Bunc)
                         O <- diag(dim(Bunc)[1])
                         colnames(O) <- rownames(O) <- rownames(Bunc)
                         O[,setdiff(c(x,m,y,covs),c(x,m,y))] <- 0
                         Iunc <- diag(dim(Bunc)[1])
                         Bstar_unc <- t(O)%*%Bunc%*%O
                         
                         ###Indirect effect matrix 
                         Indstar_unc <- solve(Iunc-Bstar_unc)-Iunc-Bstar_unc
                         
                         for(k in seq(length(umodels))){
                              xs <- x[!x %in% x[[i]]]
                              UpsES$IE$Uncond[[i]] <- Indstar_unc[xmodel.sp[[1]],xs[k]]
                              names(UpsES$IE$Uncond[[i]]) <- paste0(c(xs[k],parallel.med[[1]],y),collapse='->')
                              if(!is.null(covs)){
                                   names(UpsES$IE$Uncond[[iecount]]) <- paste0(c(names(UpsES$Uncond$Specific[iecount]),covs),collapse='|')
                              }
                         }
                         ups_unc <- Indstar_unc%*%covmatunc%*%t(Indstar_unc)
                         UpsES$Upsilon$Uncond[[i]] <- ups_unc[xmodel.sp[[1]],xmodel.sp[[1]]]
                         if(bx>1){
                              names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(paste0(x[!x %in% x[[i]]],collapse='+'),parallel.med[[1]],y),collapse='->')
                         } else names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(paste0(x[!x %in% x[[i]]],collapse='+'),parallel.med[[1]],y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(names(UpsES$Upsilon$Uncond[upscount+(i-1)*upscount]),covs),collapse='|')
                         }
                    }
                    for(i in seq(bm)){
                         O <- diag(bx+bm+length(covs)+1)
                         colnames(O) <- rownames(O) <- rownames(B)
                         O[,setdiff(c(x,m,y),c(x,m[[i]],y))] <- 0
                         Bstar_sp <- O%*%B%*%O
                         I <- diag(dim(Bstar_sp)[1])
                         Indstar_sp <- solve(I-Bstar_sp)-I-Bstar_sp
                         
                         for(k in seq(bx)){
                              UpsES$IE$Specific[[bx+k+(i-1)*bx]] <- Indstar_sp[y,x[[k]]]
                              names(UpsES$IE$Specific[[bx+k+(i-1)*bx]]) <- paste0(c(x[[k]],m[[i]],y),collapse='->')
                         }
                         if(!is.null(covs)){
                              names(UpsES$IE$Specific[[iecount]]) <- paste0(c(names(UpsES$IE$Specific[iecount]),covs),collapse='|')
                         }
                         
                         ups_sp <- Indstar_sp%*%covmat%*%t(Indstar_sp)
                         UpsES$Upsilon$Specific[[bx+i]] <- ups_sp[y,y]
                         names(UpsES$Upsilon$Specific[[bx+i]]) <- paste0(c(paste0(x,collapse='+'),m[[i]],y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$Upsilon$Specific[[bx+i]]) <- paste0(c(names(UpsES$Upsilon$Specific[upscount]),covs),collapse='|')
                         }
                         umodels <- list()
                         for(k in seq(length(models))){
                              umod <- unlist(strsplit(unlist(strsplit(models[[k]],'<-',fixed=TRUE)),'+',fixed=TRUE))
                              umod <- umod[!umod %in% m[[i]]]
                              umodels <- append(umodels,paste0(c(umod),collapse='<-'))
                         }
                         ulav.models <- list()
                         for(k in seq(length(lav.models))){
                              ulmod <- unlist(strsplit(unlist(strsplit(lav.models[[k]],'~',fixed=TRUE)),'+',fixed=TRUE))
                              ulmod <- ulmod[!ulmod %in% m[[i]]]
                              if(length(which(ulmod %in% m)>0)){
                                   if(!any(paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~')==ulav.models)){
                                        ulav.models <- append(ulav.models,paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~'))
                                   }
                              }
                         }
                         ulav.model <- paste0(rev(ulav.models),collapse='\n ')
                         out <- lavaan::sem(ulav.model,data=data)
                         ###covariance matrix from lavaan output
                         covmatunc <- lavaan::inspect(out,'cor.all')
                         Bunc <- lavaan::inspect(out,'coef')$beta
                         Bunc <- unclass(Bunc)
                         O <- diag(dim(Bunc)[1])
                         colnames(O) <- rownames(O) <- rownames(Bunc)
                         O[,setdiff(c(x,m,y,covs),c(x,m,y))] <- 0
                         Iunc <- diag(dim(Bunc)[1])
                         Bstar_unc <- t(O)%*%Bunc%*%O
                         
                         ###Indirect effect matrix 
                         Indstar_unc <- solve(Iunc-Bstar_unc)-Iunc-Bstar_unc
                         for(k in seq(bx)){
                              UpsES$IE$Uncond[[bx+k+(i-1)*bx]] <- Indstar_unc[y,x[k]]
                              names(UpsES$IE$Uncond[[bx+k+(i-1)*bx]]) <- paste0(c(x[k],paste0(m[!m %in% m[[i]]],collapse='+'),y),collapse='->')
                              if(!is.null(covs)){
                                   names(UpsES$IE$Uncond[[bx+k+(i-1)*bx]]) <- paste0(c(names(UpsES$Uncond$Specific[bx+k+(i-1)*bx]),covs),collapse='|')
                              }
                         } 
                         ups_unc <- Indstar_unc%*%covmatunc%*%t(Indstar_unc)
                         UpsES$Upsilon$Uncond[[bx+i]] <- ups_unc[y,y]
                         names(UpsES$Upsilon$Uncond[[bx+i]]) <- paste0(c(paste0(x,collapse='+'),paste0(m[!m %in% m[[i]]],collapse='+'),y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$Upsilon$Uncond[[bx+i]]) <- paste0(c(names(UpsES$Upsilon$Uncond[bx+i]),covs),collapse='|')
                         }
                    }
               } else {
                    for(i in seq(bm)){
                         O <- diag(bx+bm+length(covs)+1)
                         colnames(O) <- rownames(O) <- rownames(B)
                         O[,setdiff(c(x,m,y),c(x,m[[i]],y))] <- 0
                         Bstar_sp <- O%*%B%*%O
                         I <- diag(dim(Bstar_sp)[1])
                         Indstar_sp <- solve(I-Bstar_sp)-I-Bstar_sp
                         
                         UpsES$IE$Specific[[i]] <- Indstar_sp[y,x]
                         names(UpsES$IE$Specific[[i]]) <- paste0(c(x,m[[i]],y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$IE$Specific[[i]]) <- paste0(c(names(UpsES$IE$Specific[i]),covs),collapse='|')
                         }
                         
                         ups_sp <- Indstar_sp%*%covmat%*%t(Indstar_sp)
                         UpsES$Upsilon$Specific[[i]] <- ups_sp[y,y]
                         names(UpsES$Upsilon$Specific[[i]]) <- paste0(c(x,m[[i]],y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$Upsilon$Specific[[i]]) <- paste0(c(names(UpsES$Upsilon$Specific[i]),covs),collapse='|')
                         }
                         umodels <- list()
                         for(k in seq(length(models))){
                              umod <- unlist(strsplit(unlist(strsplit(models[[k]],'<-',fixed=TRUE)),'+',fixed=TRUE))
                              umod <- umod[!umod %in% m[[i]]]
                              umodels <- append(umodels,paste0(c(umod),collapse='<-'))
                         }
                         ulav.models <- list()
                         for(k in seq(length(lav.models))){
                              ulmod <- unlist(strsplit(unlist(strsplit(lav.models[[k]],'~',fixed=TRUE)),'+',fixed=TRUE))
                              ulmod <- ulmod[!ulmod %in% m[[i]]]
                              if(length(which(ulmod %in% m)>0)){
                                   if(!any(paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~')==ulav.models)){
                                        ulav.models <- append(ulav.models,paste0(c(ulmod[1],paste0(ulmod[ulmod!=ulmod[1]],collapse='+')),collapse='~'))
                                   }
                              }
                         }
                         ulav.model <- paste0(rev(ulav.models),collapse='\n ')
                         out <- lavaan::sem(ulav.model,data=data)
                         ###covariance matrix from lavaan output
                         covmatunc <- lavaan::inspect(out,'cor.all')
                         Bunc <- lavaan::inspect(out,'coef')$beta
                         Bunc <- unclass(Bunc)
                         O <- diag(dim(Bunc)[1])
                         colnames(O) <- rownames(O) <- rownames(Bunc)
                         O[,setdiff(c(x,m,y,covs),c(x,m,y))] <- 0
                         Iunc <- diag(dim(Bunc)[1])
                         Bstar_unc <- t(O)%*%Bunc%*%O
                         
                         ###Indirect effect matrix 
                         Indstar_unc <- solve(Iunc-Bstar_unc)-Iunc-Bstar_unc
                         for(k in seq(bx)){
                              UpsES$IE$Uncond[[i]] <- Indstar_unc[y,x]
                              names(UpsES$IE$Uncond[[i]]) <- paste0(c(x,paste0(m[!m %in% m[[i]]],collapse='+'),y),collapse='->')
                              if(!is.null(covs)){
                                   names(UpsES$IE$Uncond[[i]]) <- paste0(c(names(UpsES$Uncond$Specific[i]),covs),collapse='|')
                              }
                         } 
                         ups_unc <- Indstar_unc%*%covmatunc%*%t(Indstar_unc)
                         UpsES$Upsilon$Uncond[[i]] <- ups_unc[y,y]
                         names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(x,paste0(m[!m %in% m[[i]]],collapse='+'),y),collapse='->')
                         if(!is.null(covs)){
                              names(UpsES$Upsilon$Uncond[[i]]) <- paste0(c(names(UpsES$Upsilon$Uncond[i]),covs),collapse='|')
                         }
                    }
               }
          }
     }
     return(UpsES)
}
