# divPart development version
# includes improved performance for pairwise calculations

# Kevin Keenan 2013

# divPart, a wrapper function for the calculation of differentiation stats.
#' @export
fastDivPart <- function(infile = NULL, outfile = NULL, gp = 3, 
                        pairwise = FALSE, fst = FALSE, bs_locus = FALSE,
                        bs_pairwise = FALSE, boots = 0, plot = FALSE,
                        para = FALSE){
  
  # define arguments for testing
  #   D <- "MyData_GP2D.gen"
  #   on <- NULL
  #   fst <- T
  #   bstrps <- 10
  #   bsls <- T
  #   bspw <- T
  #   plt <- F
  #   para <- T
  #   pWise <- T
  #   gp = 2
  # define
  D <- infile
  on <- outfile
  #fst <- WC_Fst
  bstrps <- boots
  bsls <- bs_locus
  bspw <- bs_pairwise
  plt <- plot
  #para <- parallel
  pWise <- pairwise
  gp = gp
  ##############################################################################
  if(bsls==T && bstrps<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else if (bspw==T && bstrps<2){
    bs_warning <- {paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else {
    #Use pre.div to calculate the standard global and locus stats
    accDat <- pre.divLowMemory(x <- list(infile = D,
                                         gp = gp,
                                         bootstrap = FALSE,
                                         locs = TRUE,
                                         fst = fst,
                                         min = FALSE))
    # create a directory for output
    if(!is.null(on)){
      suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                             "-[diveRsity]","/",sep="")))
    }
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
    wd <- getwd()
    write_res <- is.element("xlsx", installed.packages()[, 1])
    plot_res <- is.element("sendplot", installed.packages()[, 1])
    
    para_pack_inst<-is.element(c("parallel","doParallel","foreach","iterators"),
                               installed.packages()[,1])
    
    if(plt == TRUE && is.null(on)){
      writeWarn <- paste("", "[NOTE]",
                         "Your results can't be plotted as you have not",
                         "provided an argument for 'outfile'.",
                         "Analysis completed", sep="\n")
      cat(noquote(writeWarn))
    }
    para_pack <- all(para_pack_inst)
    if(write_res == FALSE && !is.null(on)){
      Warning1<-{paste(" "," ",
                       "[NOTE]",
                       "___________________________________________________________",
                       "Please install the package 'xlsx' if you would like your", 
                       "results written to an Excel workbook.",
                       "Alternatively, your result will automatically be written",
                       "to .txt files.",
                       "___________________________________________________________",
                       "To install 'xlsx' use:",
                       "> install.packages('xlsx', dependencies=TRUE)",
                       "See:",
                       "> ?install.packages - for usage details.",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning1))
    } 
    if(plot_res==F && plt==T && !is.null(on)){
      Warning2<-{paste(" "," "," ",
                       "[NOTE]  ",
                       "___________________________________________________________",
                       "Please install the package 'sendplot' to plot your results.",
                       "Use:",
                       "> install.packages('sendplot', dependencies = TRUE)",
                       "See:",
                       "> ?install.packages - for usage details",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning2))
    }
    if(fst == TRUE){
      namer<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
               "D_Jost_est","Fst_WC","Fit_WC")
    } else {
      namer<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
               "D_Jost_est")
    }
    
    ############################################################################
    # output file multilocus stats vector 
    # pre output table for global locus stats
    
    #standard
    pre_ot1 <- cbind(accDat$locus_names, round(as.numeric(accDat$hst), 4),
                     round(as.numeric(accDat$dst), 4),
                     round(as.numeric(accDat$gst), 4),
                     round(as.numeric(accDat$gst_hedrick), 4),
                     round(as.numeric(accDat$djost), 4))
    # Add global multi locus stats to output table
    ot1 <- rbind(pre_ot1, c("Global", "", "", accDat$gst_all, 
                            accDat$gst_all_hedrick, 
                            accDat$djost_all))
    colnames(ot1) <- c("loci", "H_st", "D_st", "G_st", "G_hed_st", "D_jost")
    #Estimated
    pre_ot2 <- cbind(accDat$locus_names,
                     round(as.numeric(accDat$locus_harmonic_N),4),
                     round(as.numeric(accDat$hst_est),4),
                     round(as.numeric(accDat$dst_est),4),
                     round(as.numeric(accDat$gst_est),4),
                     round(as.numeric(accDat$gst_est_hedrick),4),
                     round(as.numeric(accDat$djost_est),4))
    
    ot2 <- rbind(pre_ot2, c("Global", "", "", "", accDat$gst_est_all, 
                            accDat$gst_est_all_hedrick, 
                            accDat$djost_est_all))
    colnames(ot2) <- c("loci", "Harmonic_N", "H_st_est", "D_st_est",
                       "G_st_est", "G_hed_st_est", "D_Jost_est")
    if(fst == TRUE){
      ot2 <- cbind(ot2, accDat$fstats[, 2:3])
    }
    if(fst == TRUE){
      plot_data321 <- c("Overall","","","",accDat$gst_est_all,
                        accDat$gst_est_all_hedrick,
                        accDat$djost_est_all,
                        as.numeric(accDat$fstats["All",2]))
      
    } else {
      plot_data321<-c("Overall","","","",accDat$gst_est_all,
                      accDat$gst_est_all_hedrick,
                      accDat$djost_est_all)
    }
    if (!is.null(on)){
      if(write_res==TRUE){
        # write data to excel
        # standard stats
        xlsx::write.xlsx(ot1,file=paste(of,"[fastDivPart].xlsx",sep=""),
                   sheetName="Standard_stats",col.names=T,
                   row.names=F,append=F)
        # Estimated stats
        xlsx::write.xlsx(ot2,file=paste(of,"[fastDivPart].xlsx",sep=""),
                   sheetName="Estimated_stats",col.names=T,
                   row.names=F,append=T)
      } else {
        # text file alternatives
        std<-file(paste(of,"Standard-stats[fastDivPart].txt",sep=""), "w")
        cat(paste(colnames(ot1),sep=""),"\n",sep="\t",file=std)
        for(i in 1:nrow(ot1)){
          cat(ot1[i,],"\n",file=std,sep="\t")
        }
        close(std)
        est<-file(paste(of,"Estimated-stats[fastDivPart].txt",sep=""),"w")
        cat(paste(colnames(ot2),sep=""),"\n",sep="\t",file=est)
        for(i in 1:nrow(ot2)){
          cat(ot2[i,],"\n",file=est,sep="\t")
        }
        close(est)
      }
    }
    ot1out<-ot1[,-1]
    ot2out<-ot2[,-1]
    
    ot1out<-matrix(as.numeric(ot1[,2:6]),ncol=5)
    rownames(ot1out)<-ot1[,1]
    colnames(ot1out)<-colnames(ot1)[-1]
    
    ot2out<-matrix(as.numeric(ot2[,-1]),ncol=(ncol(ot2)-1))
    rownames(ot2out)<-ot2[,1]
    colnames(ot2out)<-colnames(ot2)[-1]
    if (para && !para_pack){
      Warning3<-{paste(" "," ",
                       "[NOTE]",
                       "___________________________________________________________",
                       "Please make sure the packages 'parallel', 'doParallel',",
                       "'foreach' and 'iterators' are installed. These are required",
                       " to run your analysis in parallel.",
                       "Your analysis will be run sequentially!",
                       "___________________________________________________________",
                       "To install these use:",
                       "> install.packages()",
                       "See:",
                       "> ?install.packages - for usage details.",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning3))
    }
    
    ############################################################################
    ############################ Bootstrapper ##################################
    ############################################################################
    # Used only if bootstraps is greater than zero
    if(bsls == TRUE){
      
      if (para && para_pack) {
        
        if (para && para_pack) {
          #count cores
          cores <- parallel::detectCores()
          cl <- parallel::makeCluster(cores)
        }
        
        #vectorize prallele#
        gp_inls <- list(infile = D, gp = gp,
                        bootstrap = TRUE, 
                        locs = TRUE, fst = fst)
        # silence for memory efficiency
        #gp_in <- list()
        #for(i in 1:bstrps){
        #  gp_in[[i]] <- gp_inls
        #}
        
        # calculate stats from readGenepopX objects
        # export objects for parallel
        parallel::clusterExport(cl, c("gp_inls", "pre.divLowMemory"), 
                                envir = environment())
        # run parallel code
        bs_loc <- parallel::parLapply(cl, 1:bstrps, function(...){
          pre.divLowMemory(gp_inls)
        })
        # close the cluster connection
        parallel::stopCluster(cl)
        
        
        #vectorize data extraction#
        if(fst==TRUE){
          bs_glb <- do.call("rbind", lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all, 4),
              round(bs_loc[[x]]$gst_all_hedrick, 4),
              round(bs_loc[[x]]$djost_all, 4),
              round(bs_loc[[x]]$gst_est_all, 4),
              round(bs_loc[[x]]$gst_est_all_hedrick, 4),
              round(bs_loc[[x]]$djost_est_all, 4),
              as.numeric(bs_loc[[x]]$fstats["All", 2:3]))
          }))
        } else {
          bs_glb <- do.call("rbind", lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all, 4),
              round(bs_loc[[x]]$gst_all_hedrick, 4),
              round(bs_loc[[x]]$djost_all, 4),
              round(bs_loc[[x]]$gst_est_all, 4),
              round(bs_loc[[x]]$gst_est_all_hedrick, 4),
              round(bs_loc[[x]]$djost_est_all, 4))
          }))
        }
        bs_std <- lapply(1:accDat$nloci, function(x){
          do.call("rbind", lapply(1:length(bs_loc), function(y){
            c(round(bs_loc[[y]]$gst[x], 4),
              round(bs_loc[[y]]$gst_hedrick[x], 4),
              round(bs_loc[[y]]$djost[x], 4))
          }))
        })
        if(fst==TRUE){
          bs_est <- lapply(1:accDat$nloci, function(x){
            do.call("rbind", lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x], 4),
                round(bs_loc[[y]]$gst_est_hedrick[x], 4),
                round(bs_loc[[y]]$djost_est[x], 4),
                as.numeric(bs_loc[[y]]$fstats[x, 2:3]))
            }))
          })
        } else {
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4))
            }))
          })
        }
        rm(bs_loc)                  ###
        z<-gc(reset=T)                ### tidy up
        rm(z)                       ###
        
      } else {
        #vectorize non-parallel#
        
        gp_inls <- list(infile = D,
                        gp = gp,
                        bootstrap = TRUE, 
                        locs = TRUE, 
                        fst = fst)
        #gp_in<-list()
        #for(i in 1:bstrps){
        # gp_in[[i]]<-gp_inls
        #}
        # calculate stats from readGenepopX objects
        bs_loc <- lapply(1:bstrps, function(...){
          pre.divLowMemory(gp_inls)
        })
        
        
        if(fst==TRUE){
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4),
              as.numeric(bs_loc[[x]]$fstats[(accDat$nloci+1),2:3]))
          }))
        }else{
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4))
          }))
        }
        bs_std<-lapply(1:accDat$nloci, function(x){
          do.call("rbind",lapply(1:length(bs_loc), function(y){
            c(round(bs_loc[[y]]$gst[x],4),
              round(bs_loc[[y]]$gst_hedrick[x],4),
              round(bs_loc[[y]]$djost[x],4))}))
        })
        if(fst==TRUE){
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4),
                as.numeric(bs_loc[[y]]$fstats[x,2:3]))
            }))
          })
        } else {
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4))
            }))
          })
        }
        rm(bs_loc)
        z<-gc(reset=T)
        rm(z)
        
      }
      
      
      #vectorize#
      if(fst == TRUE){
        bs_res <- lapply(1:8, function(x){
          matrix(ncol = 3, nrow = (accDat$nloci+1))
        })
      } else {
        bs_res<-lapply(1:6,function(x){matrix(ncol=3, nrow=(accDat$nloci+1))})
      }
      bs_join<-cbind(bs_std, bs_est)
      bs_cis <- apply(bs_join, 1, function(x){
        res <- lapply(x, function(y){
          apply(y, 2, function(z){
            ci <- as.vector(quantile(z, probs = c(0.025, 0.975), na.rm = TRUE))
            means <- mean(z, na.rm = TRUE)
            
            return(c(means, ci))
          })
        })
        ciM <- c(res$bs_std[1,], res$bs_est[1,])
        lci <- c(res$bs_std[2,], res$bs_est[2,])
        uci <- c(res$bs_std[3,], res$bs_est[3,])
        list(mu = ciM,
             lci = lci,
             uci = uci)
      })
      mu <- t(sapply(1:length(bs_cis), function(i){
        return(bs_cis[[i]]$mu)
      }))
      lci <- t(sapply(1:length(bs_cis), function(i){
        return(bs_cis[[i]]$lci)
      }))
      uci <- t(sapply(1:length(bs_cis), function(i){
        return(bs_cis[[i]]$uci)
      }))
      # calculate ci for global
      glb_mu <- apply(bs_glb, 2, function(x){
        return(mean(x, na.rm = TRUE))
      })
      glb_lci <- apply(bs_glb, 2, function(x){
        return(quantile(x, probs = 0.025, na.rm = TRUE))
      })
      glb_uci <- apply(bs_glb, 2, function(x){
        return(quantile(x, probs = 0.975, na.rm = TRUE))
      })
      # add glb ci to mu,  uci and lci
      mu <- rbind(mu, glb_mu)
      lci <- rbind(lci, glb_lci)
      uci <- rbind(uci, glb_uci)
      #ciCalc <- function(x){
      #  res <- lapply(x, function(y){
      #    apply(y, 2, function(z){
      #      return(quantile(z, probs = c(0.025, 0.975)))
      #    })
      #  })
      #  return(res)
      #}
      #ci <- function(x){
      #  (sd(na.omit(x))/sqrt(length(na.omit(x)))) * 1.96
      #}
      #bs_cis <- t(apply(bs_join, 1, ciCalc))
      #bs_cis<-rbind(bs_cis, apply(bs_glb, 2, ci))
      if(fst==TRUE){
        for(i in 1:8){
          bs_res[[i]][,1] <- round(mu[,i], 4)
          bs_res[[i]][,2] <- round(lci[,i], 4)
          bs_res[[i]][,3] <- round(uci[,i], 4)
          bs_res[[i]][is.na(bs_res[[i]])] <- 0
        }
      } else {
        for(i in 1:6){
          bs_res[[i]][,1] <- round(mu[,i], 4)
          bs_res[[i]][,2] <- round(lci[,i], 4)
          bs_res[[i]][,3] <- round(uci[,i], 4)
          bs_res[[i]][is.na(bs_res[[i]])] <- 0
        }
      }
      
      names(bs_res) <- namer
      
      bs_res1 <- bs_res
      if(fst){
        for(i in 1:8){
          dimnames(bs_res1[[i]])<-list(c(accDat$locus_names, "global"),
                                       c("Mean","Lower_CI", "Upper_CI"))
        }
      } else {
        for(i in 1:6){
          dimnames(bs_res1[[i]])<-list(c(accDat$locus_names,"global"),
                                       c("Mean","Lower_CI","Upper_CI"))
        }
      }
      # bs results output object header
      hdr <- matrix(c("locus", "Mean", "Lower_95%CI", "Upper_95%CI"), 
                    ncol=4)
      bs_out <- matrix(rbind(hdr, c(names(bs_res)[1], "", "", ""),
                             cbind(c(accDat$locus_names, "Overall"),
                                   bs_res[[1]])), ncol = 4)
      
      if(fst){
        for(i in 2:8){
          bs_out <- matrix(rbind(bs_out, c(names(bs_res)[i], "", "", ""),
                                 cbind(c(accDat$locus_names, "global"),
                                       bs_res[[i]])), ncol = 4)
        }
      } else {
        for(i in 2:6){
          bs_out<-matrix(rbind(bs_out,c(names(bs_res)[i],"","",""),
                               cbind(c(accDat$locus_names,"Global"),
                                     bs_res[[i]])),ncol=4)
        }
      }
      if(!is.null(on)){
        if(write_res==TRUE){
          xlsx::write.xlsx(bs_out,file=paste(of,"[fastDivPart].xlsx",sep=""),
                     sheetName="Locus_bootstrap",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          bts<-file(paste(of,"Locus-bootstrap[fastDivPart].txt",sep=""), "w")
          cat(paste(colnames(bs_out),sep=""),"\n",sep="\t",file=bts)
          for(i in 1:nrow(bs_out)){
            cat(bs_out[i,],"\n",file=bts,sep="\t")
          }
          close(bts)
        }
      }
    }
    zzz<-gc()
    rm(zzz)
    if(plot_res==TRUE && plt==TRUE && bsls==TRUE){
      
      #vectorize#
      sorter<-function(x){
        z<-order(x[1:accDat$nloci,1],decreasing=F)
        #if(length(z) >= 200){
        #  z<-z[(length(z)-150):length(z)]
        #}
        return(z)
      }
      lso123<-lapply(bs_res, sorter)
      
      #
      names(lso123)<-namer
      plot.call_loci<-list()
      plot.extras_loci<-list()
      xy.labels_loci<-list()
      y.pos_loci<-list()
      x.pos_loci=1:accDat$nloci
      direct=of
      fn_pre_loci<-list()
      #Plot Gst_Nei
      plot.call_loci[[1]]=c("plot(bs_res[[4]][lso123[[4]],1],
                            ylim=c(0,(max(bs_res[[4]][,3])+
                            min(bs_res[[4]][,3]))),xaxt='n',
                            ylab=names(bs_res)[4],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")
      
      plot.extras_loci[[1]]=c("points(bs_res[[4]][lso123[[4]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[4]][lso123[[4]],2],
                              1:accDat$nloci,bs_res[[4]][lso123[[4]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[4]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
      
      xy.labels_loci[[1]]=data.frame(Locus_name=accDat$locus_names[lso123[[4]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[4]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[4]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[4]],1],4))
      
      y.pos_loci[[1]]=bs_res[[4]][lso123[[4]],1]
      fn_pre_loci[[1]]<-names(bs_res)[4]
      
      
      
      # Plot Gst_Hedrick
      plot.call_loci[[2]]=c("plot(bs_res[[5]][lso123[[5]],1],
                            ylim=c(0,1),xaxt='n',ylab=names(bs_res)[5],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")
      
      plot.extras_loci[[2]]=c("points(bs_res[[5]][lso123[[5]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[5]][lso123[[5]],2],
                              1:accDat$nloci,bs_res[[5]][lso123[[5]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[5]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
      
      xy.labels_loci[[2]]=data.frame(Locus_name=accDat$locus_names[lso123[[5]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[5]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[5]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[5]],1],4))
      
      y.pos_loci[[2]]=bs_res[[5]][lso123[[5]],1]
      fn_pre_loci[[2]]<-names(bs_res)[5]
      
      
      # Plot D_jost
      plot.call_loci[[3]]=c("plot(bs_res[[6]][lso123[[6]],1],
                            ylim=c(0,1),xaxt='n',ylab=names(bs_res)[6],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")
      
      plot.extras_loci[[3]]=c("points(bs_res[[6]][lso123[[6]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[6]][lso123[[6]],2],
                              1:accDat$nloci,bs_res[[6]][lso123[[6]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[6]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
      
      xy.labels_loci[[3]]=data.frame(Locus_name=accDat$locus_names[lso123[[6]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[6]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[6]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[6]],1],4))
      
      y.pos_loci[[3]]=bs_res[[6]][lso123[[6]],1]
      fn_pre_loci[[3]]<-names(bs_res)[6]
      
      #plot(Fst)
      if(fst==TRUE){
        plot.call_loci[[4]]=c("plot(bs_res[[8]][lso123[[8]],1],
                              ylim=c(0,(max(bs_res[[8]][,3])+
                              min(bs_res[[8]][,3]))),xaxt='n',
                              ylab=names(bs_res)[8],type='n',
                              xlab='Loci \n (Hover over a point to see locus data)',
                              cex.lab=1.5,cex.axis=1.3,las=1)")
        
        plot.extras_loci[[4]]=c("points(bs_res[[8]][lso123[[8]],1],
                                pch=15,col='black',cex=1);
                                arrows(1:accDat$nloci,bs_res[[8]][lso123[[8]],2],
                                1:accDat$nloci,bs_res[[8]][lso123[[8]],3],code=3,
                                angle=90,length=0.05,lwd=0.1);
                                abline(h=c(0,bs_res[[8]][(accDat$nloci+1),2]),
                                lwd=1,lty=c(1,2),col=c('black','red'))")
        
        xy.labels_loci[[4]]=data.frame(Locus_name=accDat$locus_names[lso123[[8]]],
                                       Gst_Nei=round(bs_res[[4]][lso123[[8]],1],4),
                                       Gst_Hedrick=round(bs_res[[5]][lso123[[8]],1],4),
                                       D_jost=round(bs_res[[6]][lso123[[8]],1],4),
                                       Fst_WC=round(bs_res[[8]][lso123[[8]],1],4))
        
        y.pos_loci[[4]]=bs_res[[8]][lso123[[8]],1]
        fn_pre_loci[[4]]<-names(bs_res)[8]
      }
    }
    ############################################################################
    ################################## Pairwise ################################
    ############################################################################
    # population pair combinations
    
    # define new functions
    ############################################################################
    ############################################################################
    # pwCalc
    ############################################################################
    # New optimised function for the calculation of pairwise statistics
    # Returns a 3D array where each 'slot' represents the pairwise matrix
    # for Gst_est, G_st_est_hed and D_jost_est respectively
    
    # Kevin Keenan
    # 2013
    
    pwCalc <- function(infile, fst,  bs = FALSE){
      
      
      #   # uncomment for testing
      #   infile <- "pw_test.txt"
      #   source("readGenepopX.R")
      #   # read pwBasicCalc function
      #   source("pwBasicCalc.R")
      # define baseline info
      dat <- readGenepopX(list(infile = infile,
                               bootstrap = bs))
      if(fst){
        # calculate all fst
        fstat <- pwFstWC(dat)
        # extract locus theta and variance components
        locTheta <- lapply(fstat, "[[", 1)
        # sum res
        aLoc <- Reduce(`+`, lapply(fstat, "[[", 2))
        bLoc <- Reduce(`+`, lapply(fstat, "[[", 3))
        cLoc <- Reduce(`+`, lapply(fstat, "[[", 4))
        # calculate pw Fst across loci
        pwTheta <- aLoc/(aLoc+bLoc+cLoc)
        # clean up
        rm(aLoc, bLoc, cLoc, fstat)
        z <- gc()
        rm(z)
      }
      # extract allele frequencies
      af <- dat$allele_freq
      # extract harmonic mean sample sizes
      
      # make space in RAM
      dat$allele_freq <- NULL
      z <- gc()
      rm(z)
      # extract npops and nloci
      npops <- dat$npops
      nloci <- dat$nloci
      # define pairwise relationships
      pw <- combn(dat$npops, 2)
      # generate pairwise locus harmonic mean sample sizes
      indtyp <- dat$indtyp
      pwHarm <- lapply(indtyp, pwHarmonic, pw = pw)
      
      
      # calculate pairwise ht and hs
      hths <- mapply(pwBasicCalc, af, pwHarm,
                     MoreArgs = list(pw = pw, npops = dat$npops),
                     SIMPLIFY = FALSE)
      # seperate ht and hs
      #   htLoc <- lapply(hths, "[[", 1)
      #   hsLoc <- lapply(hths, "[[", 2)
      # seperate ht_est and hs_est
      hsEstLoc <- lapply(hths, "[[", 1)
      htEstLoc <- lapply(hths, "[[", 2)
      
      # clean up
      rm(hths)
      z <- gc()
      rm(z)
      
      # Calculate locus stats
      # Standard locus stats
      # locus Gst
      #   gstLoc <- mapply(FUN = gstCalc, ht = htLoc, hs = hsLoc, 
      #                    SIMPLIFY = FALSE)
      #   # locus G'st
      #   gstHedLoc <- mapply(FUN = gstHedCalc, ht = htLoc, hs = hsLoc,
      #                       SIMPLIFY = FALSE)
      #   # locus D_jost
      #   dLoc <- mapply(FUN = djostCalc, ht = htLoc, hs = hsLoc,
      #                  SIMPLIFY = FALSE)
      
      # Estimated locus stats
      # locus Gst_est
      gstLocEst <- mapply(FUN = gstCalc, ht = htEstLoc, 
                          hs = hsEstLoc, 
                          SIMPLIFY = FALSE)
      # locus G'st_est
      gstHedLocEst <- mapply(FUN = gstHedCalc, ht = htEstLoc, 
                             hs = hsEstLoc,
                             SIMPLIFY = FALSE)
      # locus D_jost_est
      dLocEst <- mapply(FUN = djostCalc, ht = htEstLoc, 
                        hs = hsEstLoc,
                        SIMPLIFY = FALSE)
      
      #   # calculate mean ht and hs
      #   htMean <- Reduce(`+`, htLoc)/nloci
      #   hsMean <- Reduce(`+`, hsLoc)/nloci
      # calculate mean ht_est and hs_est
      htEstMean <- Reduce(`+`, htEstLoc)/nloci
      hsEstMean <- Reduce(`+`, hsEstLoc)/nloci
      
      # calculate standard stats (uncomment for loc stats)
      
      #   # overall dst
      #   dstAll <- htMean - hsMean
      #   # overall gst (Nei 1973)
      #   gstAll <- (dstAll)/htMean
      #   # overall max gst (Hedricks 2005)
      #   gstAllMax <- ((2 - 1)*(1 - hsMean)) / ((2 - 1) + hsMean)
      #   # overall Hedricks' Gst
      #   gstAllHedrick <- gstAll/gstAllMax
      #   # Overall D_jost (Jost 2008)
      #   djostAll <- (dstAll/(1-hsMean))*(2/(2-1))
      
      # Calculate estimated stats
      
      # Overall estimated dst
      dstEstAll <- htEstMean - hsEstMean
      # Overall estimated Gst (Nei & Chesser, 1983)
      gstEstAll <- dstEstAll/htEstMean
      # Overall estimated max Gst (Hedricks 2005)
      gstEstAllMax <- ((2-1)*(1-hsEstMean))/(2-1+hsEstMean)
      # Overall estimated Hedricks' Gst
      gstEstAllHed <- gstEstAll/gstEstAllMax
      # Overall estimated D_Jost (Chao et al., 2008)
      if(nloci == 1){
        djostEstAll <- (2/(2-1))*((dstEstAll)/(1 - hsEstMean))
      } else {
        dLocEstMn <- Reduce(`+`, dLocEst)/nloci
        # calculate variance (convert dLocEst to an array)
        dLocEst1 <- array(unlist(dLocEst), 
                          dim = c(nrow(dLocEst[[1]]), 
                                  ncol(dLocEst[[1]]), 
                                  length(dLocEst)))
        dLocEstVar <- apply(dLocEst1, c(1,2), var)
        djostEstAll <- 1/((1/dLocEstMn)+((dLocEstVar*((1/dLocEstMn)^3))))
        # tidy up
        rm(dLocEstMn, dLocEstVar)
        z <- gc()
        rm(z)
      }
      
      # define a function to arrange locus stats into arrays
      #   arrDef <- function(x){
      #     return(array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
      #   }
      if(fst){
        resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll, pwTheta),
                        dim = c(nrow(gstEstAll),
                                ncol(gstEstAll),
                                4))
        lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 4))
        lstats[,,,1] <- unlist(gstLocEst)
        lstats[,,,2] <- unlist(gstHedLocEst)
        lstats[,,,3] <- unlist(dLocEst)
        lstats[,,,4] <- unlist(locTheta)
        #lstats <- array(list(gstLocEst, gstHedLocEst, dLocEst, locTheta)
        #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst, locTheta,
        #                 SIMPLIFY=FALSE)
      } else {
        resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll),
                        dim = c(nrow(gstEstAll),
                                ncol(gstEstAll), 3))
        lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 3))
        lstats[,,,1] <- unlist(gstLocEst)
        lstats[,,,2] <- unlist(gstHedLocEst)
        lstats[,,,3] <- unlist(dLocEst)
        #lstats <- list(gstLocEst, gstHedLocEst, dLocEst)
        #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst,
        #                 SIMPLIFY = FALSE)
      }
      # arrange loci into arrays
      #locOut <- lapply(lstats, arrDef)
      
      
      list(resArr = resArr,
           locOut = lstats)
    }
    ############################################################################
    # END - pwDivCalc
    ############################################################################
    # Calculate Weir & Cockerham's F-statistics (optimised)
    ############################################################################
    # pwFstWC: a function co calculate weir and cockerhams fis, fit, and fst
    ############################################################################
    pwFstWC <- function(rdat){
      #   rdat <- diveRsity::readGenepop("KK_test1v2.gen")
      pw <- combn(rdat$npops, 2)
      #   # account for loci with missing info for pops
      #   pwBadData <- function(indtyp, pw){
      #     out <- sapply(1:ncol(pw), function(i){
      #       is.element(0, indtyp[pw[,i]])
      #     })
      #   }
      #   badDat <- sapply(rdat$indtyp, pwBadData, pw = pw)
      #   if(any(badDat)){
      #     bd <- TRUE
      #   }
      #   # determine the number of loci per pw comparison
      #   nlocPw <- apply(badDat, 1, function(x){
      #     if(sum(x) > 0){
      #       nl <- rdat$nloci - sum(x)
      #     } else {
      #       nl <- rdat$nloci
      #     }
      #   })
      #   # define all good data
      #   gdDat <- lapply(1:nrow(badDat), function(i){
      #     which(!badDat[i,])
      #   })
      #   badDat <- lapply(1:nrow(badDat), function(i){
      #     which(badDat[i,])
      #   })
      # get all genotypes for each pw comparison
      allGenot <- apply(pw, 2, function(x){
        list(rdat$pop_list[[x[1]]], 
             rdat$pop_list[[x[2]]])
      })
      #   # filter bad data
      #   if(any(nlocPw != rdat$nloci)){
      #     idx <- which(nlocPw != rdat$nloci)
      #     for(i in idx){
      #       allGenot[[i]][[1]] <- allGenot[[i]][[1]][, gdDat[[i]]]
      #       allGenot[[i]][[2]] <- allGenot[[i]][[2]][, gdDat[[i]]]
      #     }
      #   }
      # unlist pw genotype data
      allGenot <- lapply(allGenot, function(x){
        return(do.call("rbind", x))
      })
      # identify unique genotypes
      #   genot <- lapply(allGenot, function(x){
      #     return(apply(x, 2, function(y){
      #       unique(na.omit(y))
      #     }))
      #   })
      # count number of genotypes per pw per loc
      
      genoCount <- lapply(allGenot, function(x){
        if(NCOL(x) == 1){
          return(list(table(x)))
        } else {
          lapply(1:ncol(x), function(i) table(x[,i]))
        }
      })
      
      
      #   genoCount <- lapply(allGenot, function(x){
      #     lapply(split(x,seq(NCOL(x))),table) # accounts for single loci
      #     #apply(x, 2, table)
      #   })
      
      
      # function to count heterozygotes
      htCount <- function(x){
        nms <- names(x)
        ncharGeno <- nchar(nms[1])
        alls <- cbind(substr(nms, 1, (ncharGeno/2)),
                      substr(nms, ((ncharGeno/2) + 1), ncharGeno))
        unqAlls <- unique(as.vector(alls))
        hetCounts <- sapply(unqAlls, function(a){
          idx <- which(rowSums(alls == a) == 1)
          return(sum(x[idx]))
        })
        # needed for no het count
        if(length(hetCounts) == 0L){
          return(NA)
        } else {
          return(hetCounts) 
        }
      }
      # hSum is the total observed hets per allele
      hSum <- lapply(genoCount, function(x){
        out <- lapply(x, htCount)
      })
      
      #   if(bd){
      #     # insert na for missing loci
      #     hSum <- lapply(seq_along(badDat), function(i){
      #       naPos <- badDat[[i]]
      #       idx <- c(seq_along(hSum[[i]]), (naPos - 0.5))
      #       return(c(hSum[[i]], rep(NA, length(naPos)))[order(idx)])
      #     }) 
      #   }
      # convert to locus orientated hSum
      hSum <- lapply(seq_along(hSum[[1]]), function(i){
        lapply(hSum, "[[", i)
      })
      
      # total ind typed per loc per pw
      indTypTot <- lapply(rdat$indtyp, function(x){
        return(apply(pw, 2, function(y){
          sum(x[y], na.rm = TRUE)
        }))
      })
      # nBar is the mean number of inds per pop
      nBar <- lapply(indTypTot, `/`, 2)
      
      # hbar per pw per loc
      hBar <- lapply(seq_along(hSum), function(i){
        divd <- indTypTot[[i]]
        return(mapply(`/`, hSum[[i]], divd, SIMPLIFY = FALSE))
      })
      
      # p per loc per pw
      pCalc <- function(x, y, pw){
        out <- lapply(seq_along(pw[1,]), function(i){
          return(cbind((x[,pw[1,i]]*(2*y[pw[1,i]])),
                       (x[,pw[2,i]]*(2*y[pw[2,i]]))))
        })
        return(out)
      }
      p <- mapply(FUN = pCalc, x = rdat$allele_freq, 
                  y = rdat$indtyp, 
                  MoreArgs = list(pw = pw), 
                  SIMPLIFY = FALSE)
      
      #   # convert p elements into array structure
      #   pArr <- lapply(p, function(x){
      #     d3 <- length(x)
      #     d2 <- 2
      #     d1 <- nrow(x[[1]])
      #     return(array(unlist(x), dim = c(d1, d2, d3)))
      #   })
      
      fstatCal <- function(indT, indtyp, hBar, nBar, p, pw, npops){
        #                 indT=indTypTot[[17]]
        #                 indtyp=rdat$indtyp[[17]]
        #                 hBar <- hBar[[17]]
        #                 nBar <- nBar[[17]]
        #                 p <- p[[17]]
        #                 pw <- pw
        #                 npops <- rdat$npops
        indLocPwSqSum <- sapply(seq_along(pw[1,]), function(i){
          return(sum(indtyp[pw[,i]]^2))
        })
        indtypPw <- lapply(1:ncol(pw), function(idx){
          return(indtyp[pw[,idx]])
        })
        nC <- indT - (indLocPwSqSum/indT)
        ptildCalc <- function(x,y){ 
          return(cbind((x[,1]/(2*y[1])),
                       (x[,2]/(2*y[2]))))
        }
        pTild <- mapply(FUN = ptildCalc, x = p, y = indtypPw,
                        SIMPLIFY = FALSE)
        pBar <- lapply(seq_along(p), function(i){
          return(rowSums((p[[i]])/(2*indT[i])))
        })
        s2 <- lapply(seq_along(pBar), function(i){
          pp <- (pTild[[i]]-pBar[[i]])^2
          pp <- cbind((pp[,1]*indtypPw[[i]][1]),
                      (pp[,2]*indtypPw[[i]][2]))
          pp <- rowSums(pp)
          return((pp/(1*nBar[i])))
        })
        A <- lapply(seq_along(pBar), function(i){
          return(pBar[[i]]*(1-pBar[[i]])-(1)*s2[[i]]/2)
        })
        # fix hBar for unequal lengths
        idx <- lapply(seq_along(A), function(i){
          out <- match(names(A[[i]]), names(hBar[[i]]))
          return(which(!is.na(out)))
        })
        A <- lapply(seq_along(A), function(i){
          return(A[[i]][idx[[i]]])
        })
        s2 <- lapply(seq_along(s2), function(i){
          return(s2[[i]][idx[[i]]])
        })
        a <- lapply(seq_along(s2), function(i){
          return(nBar[[i]]*(s2[[i]]-(A[[i]]-(hBar[[i]]/4))/(nBar[[i]]-1))/nC[[i]])
        })
        #     a <- lapply(seq_along(s2), function(i){
        #       return(a[[i]][idx[[i]]])
        #     })
        b <- lapply(seq_along(A), function(i){
          return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-((2*nBar[[i]]-1)/(4*nBar[[i]]))*hBar[[i]]))
          #return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-(2*(nBar[[i]]-1))*hBar[[i]]/(4*nBar[[i]])))
        })
        #     b <- lapply(seq_along(A), function(i){
        #       return(b[[i]][idx[[i]]])
        #     })
        cdat <- lapply(seq_along(A), function(i){
          return(hBar[[i]]/2)
        })
        #     cdat <- lapply(seq_along(A), function(i){
        #       return(cdat[[i]][idx[[i]]])
        #     })
        A <- sapply(A, sum)
        a <- sapply(a, sum)
        b <- sapply(b, sum)
        cdat <- sapply(cdat, sum)
        cdat[is.na(cdat)] <- 0
        theta <- a/(a+b+cdat)
        theta[is.nan(theta)] <- NA
        pwMat <- matrix(ncol = npops, nrow = npops)
        aMat <- matrix(ncol = npops, nrow = npops)
        bMat <- matrix(ncol = npops, nrow = npops)
        cMat <- matrix(ncol = npops, nrow = npops)
        for(i in 1:ncol(pw)){
          pwMat[pw[2,i], pw[1,i]] <- theta[i]
          aMat[pw[2,i], pw[1,i]] <- a[i]
          bMat[pw[2,i], pw[1,i]] <- b[i]
          cMat[pw[2,i], pw[1,i]] <- cdat[i]
        }
        #pwMat[is.nan(pwMat)] <- 0
        aMat[is.nan(aMat)] <- 0
        cMat[is.nan(bMat)] <- 0
        bMat[is.nan(bMat)] <- 0
        
        list(pwMat, aMat, bMat, cMat)
      }
      
      # run fstatCal for each locus
      pwLoc <- mapply(FUN = fstatCal, indT = indTypTot,
                      indtyp = rdat$indtyp, hBar = hBar,
                      nBar = nBar, p = p, 
                      MoreArgs = list(pw = pw, npops = rdat$npops),
                      SIMPLIFY = FALSE)
      return(pwLoc)
    }
    ############################################################################
    # END - pwDivCalc
    ############################################################################
    # pwBasicCalc: a small function for calculating pairwise ht and hs 
    ############################################################################
    pwBasicCalc <- function(af, sHarm, pw, npops){
      ht <- matrix(ncol = npops, nrow = npops)
      hs <- matrix(ncol = npops, nrow = npops)
      htEst <- matrix(ncol = npops, nrow = npops)
      hsEst <- matrix(ncol = npops, nrow = npops)
      for(i in 1:ncol(pw)){
        id1 <- pw[1,i]
        id2 <- pw[2,i]
        # locus ht
        ht[id2, id1] <- 1 - sum(((af[,id1] + af[,id2])/2)^2)
        # locus hs
        hs[id2, id1] <- 1 - sum((af[,id1]^2 + af[,id2]^2)/2)
        # locus hs_est
        hsEst[id2, id1] <- hs[id2, id1]*((2*sHarm[id2,id1])/(2*sHarm[id2,id1]-1))
        # locus ht_est
        htEst[id2, id1] <- ht[id2, id1] + (hsEst[id2, id1]/(4*sHarm[id2, id1]))
      }
      #   ht[is.nan(ht)] <- 0
      #   hs[is.nan(hs)] <- 0
      htEst[is.nan(htEst)] <- 0
      hsEst[is.nan(hsEst)] <- 0
      list(hsEst = hsEst,
           htEst = htEst)
    }
    ############################################################################
    # END - pwBasicCalc
    ############################################################################
    
    # define locus stat calculators
    gstCalc <- function(ht, hs){
      return((ht - hs)/ht)
    }
    
    gstHedCalc <- function(ht, hs){
      gstMax <- ((2-1)*(1-hs))/(2-1+hs)
      return(((ht-hs)/ht)/gstMax)
    }
    
    djostCalc <- function(ht, hs){
      return((2/1)*((ht-hs)/(1-hs)))
    }
    
    # calculate pairwise locus harmonic mean
    pwHarmonic <- function(lss, pw){
      np <- length(lss)
      lhrm <- matrix(ncol = np, nrow = np)
      pwSS <- cbind(lss[pw[1,]], lss[pw[2,]])
      lhrmEle <- (0.5 * ((pwSS[,1]^-1) + (pwSS[,2]^-1)))^-1
      for(i in 1:ncol(pw)){
        idx1 <- pw[1,i]
        idx2 <- pw[2,i]
        lhrm[idx2, idx1] <- lhrmEle[i]
      }
      return(lhrm)
    }
    ############################################################################
    # pwDivCalc: a small function for calculating pairwise ht and hs 
    ############################################################################
    pwDivCalc <- function(x, pw, npops){
      ht <- matrix(ncol = npops, nrow = npops)
      hs <- matrix(ncol = npops, nrow = npops)
      for(i in 1:ncol(pw)){
        gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
        f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
        ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
        ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
        hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
        hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
      }
      ht[is.nan(ht)] <- 0
      hs[is.nan(hs)] <- 0
      list(ht = ht, 
           hs = hs)
    }
    ############################################################################
    # END - pwDivCalc
    ############################################################################
    
    ############################################################################
    ############################################################################
    # working well 24/10/13
    if(pWise || bspw){
      # get pw names
      pw <- combn(accDat$npops, 2)
      popNms <- accDat$pop_names
      # for pw bootstrap table
      pw_nms <- paste(popNms[pw[1,]], popNms[pw[2,]], sep = " vs. ")
      
      pwStats <- pwCalc(D, fst, bs = FALSE)
      # extract stats
      gstPW <- pwStats$resArr[,,1]
      gstHPW <- pwStats$resArr[,,2]
      dPW <- pwStats$resArr[,,3]
      if(fst){
        thetaPW <- pwStats$resArr[,,4]
      }
      # clean up
      locstats <- pwStats$locOut
      rm(pwStats)
      z <- gc()
      rm(z)
      spc1 <- rep("", ncol(gstPW))
      if(fst){
        statNms <- c("Gst_est", "G'st_est", "Djost_est", "Fst_WC")
        outobj <- rbind(c(statNms[1], spc1), 
                        c("", popNms),
                        cbind(popNms, round(gstPW, 4)),
                        c(statNms[2], spc1),
                        c("", popNms),
                        cbind(popNms, round(gstHPW, 4)), 
                        c(statNms[3], spc1),
                        c("", popNms),
                        cbind(popNms, round(dPW, 4)), 
                        c(statNms[4], spc1),
                        c("", popNms),
                        cbind(popNms, round(thetaPW, 4)))
        outobj[is.na(outobj)] <- ""
        pwMatListOut <- list(gstPW, gstHPW, dPW, thetaPW)
        # add names to pwMatListOut
        names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
        # tidy up
        rm(gstPW, gstHPW, dPW, thetaPW)
        z <- gc()
        rm(z)
      } else {
        statNms <- c("Gst_est", "G'st_est", "Djost_est")
        outobj <- rbind(c(statNms[1], spc1), 
                        c("", popNms),
                        cbind(popNms, round(gstPW, 4)),
                        c(statNms[2], spc1),
                        c("", popNms),
                        cbind(popNms, round(gstHPW, 4)), 
                        c(statNms[3], spc1),
                        c("", popNms),
                        cbind(popNms, round(dPW, 4)))
        outobj[is.na(outobj)] <- ""
        pwMatListOut <- list(gstPW, gstHPW, dPW)
        # add names to pwMatListOut
        names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst")
        # tidy up
        rm(gstPW, gstHPW, dPW)
        z <- gc()
        rm(z)
      }
      if(!is.null(on)){
        if(write_res == TRUE){
          # write data to excel
          # Load dependencies
          # pw stats
          xlsx::write.xlsx(outobj, file = paste(of, "[fastDivPart].xlsx", 
                                                sep=""),
                     sheetName = "Pairwise-stats", col.names = FALSE,
                     row.names = FALSE, append = TRUE)
        } else {
          # text file alternatives
          pw_outer <- file(paste(of, "Pairwise-stats[fastDivPart].txt", 
                                 sep=""), "w")
          for(i in 1:nrow(outobj)){
            cat(outobj[i,], "\n", file = pw_outer, sep = "\t")
          }
          close(std)
        }
      }
      for(i in 1:length(pwMatListOut)){
        dimnames(pwMatListOut[[i]]) <- list(popNms, popNms)
      }
      # convert locstats in to list format
      locstats <- lapply(apply(locstats, 4, list), function(x){
        lapply(apply(x[[1]], 3, list), function(y){
          out <- y[[1]]
          dimnames(out) <- list(popNms, popNms)
          return(out)
        })
      })
      # add names etc
      if(fst){
        # prepare locstats for output
        names(locstats) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
      } else {
        names(locstats) <- c("gstEst", "gstEstHed", "djostEst")
      }
      for(i in 1:length(locstats)){
        names(locstats[[i]]) <- accDat$locus_names
      }
      
    }
    
    #Bootstrap
    if(bspw == TRUE){
      if (para && para_pack) {
        
        cl <- parallel::makeCluster(detectCores())
        parallel::clusterExport(cl, c("pwCalc", "fst", "D", "readGenepopX",
                                      "fileReader", "pwFstWC", "pwHarmonic",
                                      "pwBasicCalc", "djostCalc", "gstCalc",
                                      "gstHedCalc"), 
                                envir = environment())
        pwBsStat <- parallel::parLapply(cl, 1:bstrps, function(...){
          return(pwCalc(infile = D, fst, bs = TRUE))
        })
        parallel::stopCluster(cl)
        
      } else {
        pwBsStat <- lapply(1:bstrps, function(...){
          return(pwCalc(D, fst, bs = TRUE))
        })
      }
      
      
      # seperate each stat
      
      gstEst <- lapply(pwBsStat, function(x){
        x$resArr[,,1]
      })
      
      gstEstHed <- lapply(pwBsStat, function(x){
        x$resArr[,,2]
      })
      
      dEst <- lapply(pwBsStat, function(x){
        x$resArr[,,3]
      })
      
      if(fst){
        theta <- lapply(pwBsStat, function(x){
          x$resArr[,,4]
        })
      }
      pwBsLoc <- lapply(pwBsStat, "[[", 2)
      # tidy up
      rm(pwBsStat)
      z <- gc()
      rm(z)
      
      # convert bs lists to arrays for calculations
      if(fst){
        stats <- list(gstEst = array(unlist(gstEst),
                                     dim = c(nrow(gstEst[[1]]),
                                             nrow(gstEst[[1]]),
                                             bstrps)),
                      gstEstHed = array(unlist(gstEstHed),
                                        dim = c(nrow(gstEstHed[[1]]),
                                                nrow(gstEstHed[[1]]),
                                                bstrps)),
                      dEst = array(unlist(dEst),
                                   dim = c(nrow(dEst[[1]]),
                                           nrow(dEst[[1]]),
                                           bstrps)),
                      theta = array(unlist(theta),
                                    dim = c(nrow(theta[[1]]),
                                            nrow(theta[[1]]),
                                            bstrps)))
        
      } else {
        stats <- list(gstEst = array(unlist(gstEst),
                                     dim = c(nrow(gstEst[[1]]),
                                             nrow(gstEst[[1]]),
                                             bstrps)),
                      gstEstHed = array(unlist(gstEstHed),
                                        dim = c(nrow(gstEstHed[[1]]),
                                                nrow(gstEstHed[[1]]),
                                                bstrps)),
                      dEst = array(unlist(dEst),
                                   dim = c(nrow(dEst[[1]]),
                                           nrow(dEst[[1]]),
                                           bstrps)))
      }
      # tidy up
      if(fst){
        rm(dEst, gstEst, gstEstHed, theta)
        z <- gc()
        rm(z) 
      } else {
        # tidy up
        z <- gc()
        rm(z) 
      }
      # convert locus stats into arrays for CI calculations
      npops <- accDat$npops
      nloci <- accDat$nloci
      if(fst){
        locStats <- list(
          # nei's Gst
          gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,1])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Hedrick's Gst
          gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,2])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Jost's D
          dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,3])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Weir & Cockerham's Fst
          thetaLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,4])
          })), dim = c(npops, npops, nloci, bstrps))
        )
      } else {
        locStats <- list(
          # nei's Gst
          gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,1])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Hedrick's Gst
          gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,2])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Jost's D
          dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,3])
          })), dim = c(npops, npops, nloci, bstrps))
        )
      }
      # convert locus bs stats into seperate loci
      locStats <- lapply(locStats, function(x){
        lapply(apply(x, 3, list), function(y){
          return(y[[1]])
        })
      })
      
      # calculate bias corrected CI
      
      biasCor <- function(param, bs_param){
        mnBS <- apply(bs_param, c(1,2), mean, na.rm = TRUE)
        mnBS[is.nan(mnBS)] <- NA
        mnBS <- mnBS - param
        bs_param <- sweep(bs_param, c(1:2), mnBS, "-")
        return(bs_param)
      }
      biasFix <- lapply(1:length(locstats), function(i){
        mapply(FUN = biasCor, param = locstats[[i]], bs_param = locStats[[i]],
               SIMPLIFY = FALSE)
      })
      # works well
      if (para && para_pack) {
        
        cl <- parallel::makeCluster(detectCores())
        bcLocLCI <- parallel::parLapply(cl, biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
          })
          return(loc)
        })
        bcLocUCI <- parallel::parLapply(cl, biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
          })
          return(loc)
        })
        parallel::stopCluster(cl)
      } else {
        bcLocLCI <- lapply(biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
          })
          return(loc)
        })
        bcLocUCI <- lapply(biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
          })
          return(loc)
        })
      }
      # organise data into output format
      # define function
      dfSort <- function(act, Low, High, pw_nms){
        df <- data.frame(actual = as.vector(act[lower.tri(act)]),
                         lower = as.vector(Low[lower.tri(Low)]),
                         upper = as.vector(High[lower.tri(High)]))
        rownames(df) <- pw_nms
        df[is.nan(as.matrix(df))] <- NA
        return(df)
      }
      pwLocOutput <- lapply(1:length(locstats), function(i){
        mapply(FUN = dfSort, act = locstats[[i]], Low = bcLocLCI[[i]],
               High = bcLocUCI[[i]], MoreArgs = list(pw_nms = pw_nms),
               SIMPLIFY = FALSE)
      })
      if(fst){
        names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
      } else {
        names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst")
      }
      for(i in 1:length(pwLocOutput)){
        names(pwLocOutput[[i]]) <- accDat$locus_names
      }
      
      
      
      # uncomment for standard CIs
      #       # calculate the CIs per locus
      #       if (para && para_pack) {
      #         library(parallel)
      #         cl <- makeCluster(detectCores())
      #         locLCI <- parLapply(cl, locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #         # upper CIs
      #         locUCI <- parLapply(cl, locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #       } else {
      #         locLCI <- lapply(locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #         # upper CIs
      #         locUCI <- lapply(locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #       }
      
      #         
      #       # organise locus stats into dataframe
      #       # pw sorter
      #       pwSorter <- function(x, pw){
      #         return(data.frame(actual = sapply(1:ncol(pw), function(i){
      #           x[[1]][pw[2,i], pw[1,i]]
      #           }), 
      #           Lower = sapply(1:ncol(pw), function(i){
      #             x[[2]][pw[2,i], pw[1,i]]
      #           }),
      #           Upper = sapply(1:ncol(pw), function(i){
      #             x[[3]][pw[2,i], pw[1,i]]
      #           })
      #         ))
      #       }
      #       locOutCol <- lapply(locOut, function(x){
      #         lapply(x, pwSorter, pw = pw)
      #       })
      
      
      
      
      
      
      # organise data
      # calculate bias for cis
      biasCalc <- function(param, bs_param, pw){
        #bias <- param
        for(i in 1:ncol(pw)){
          dat <- bs_param[pw[2,i], pw[1,i], ]
          t0 <- param[pw[2,i], pw[1,i]]
          mnBS <- mean(dat , na.rm = TRUE) - t0
          bs_param[pw[2,i], pw[1,i], ] <- bs_param[pw[2,i], pw[1,i], ] - mnBS
        }
        return(bs_param)
      }
      
      # try adjusting bootstrapped estimate using bias
      
      bcStats <- mapply(biasCalc, param = pwMatListOut, bs_param = stats, 
                        MoreArgs = list(pw = pw), SIMPLIFY = FALSE)
      
      # calculate the upper and lower 95% ci
      lowCI <- lapply(stats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.025, na.rm = TRUE))
      })
      
      # bias corrected
      bcLowCI <- lapply(bcStats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.025, na.rm = TRUE))
      })
      
      
      upCI <- lapply(stats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.975, na.rm = TRUE))
      })
      
      # bias corrected
      bcHighCI <- lapply(bcStats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.975, na.rm = TRUE))
      })
      
      
      statMean <- lapply(stats, function(x){
        return(apply(x, c(1,2), mean, na.rm = TRUE))
      })
      
      # bias corrected
      bcStatMean <- lapply(bcStats, function(x){
        return(apply(x, c(1,2), mean, na.rm = TRUE))
      }) 
      
      # tidy up
      rm(stats)
      z <- gc()
      rm(z)
      
      # organize ci and mean into output structure
      pw <- combn(ncol(lowCI[[1]]), 2)
      outOrg <- function(t0 ,t1 , t2, l1, l2, u1, u2, pw, pwNms){
        out <- matrix(ncol = 7, nrow = ncol(pw))
        colnames(out) <- c("actual", "mean", "BC_mean", "Lower_95%CI", 
                           "Upper_95%CI", "BC_Lower_95%CI", "BC_Upper_95%CI")
        rownames(out) <- pwNms
        for(i in 1:ncol(pw)){
          idx <- as.vector(rev(pw[,i]))
          out[i,] <- c(t0[idx[1], idx[2]], t1[idx[1], idx[2]],
                       t2[idx[1], idx[2]], l1[idx[1], idx[2]],
                       l2[idx[1], idx[2]], u1[idx[1], idx[2]],
                       u2[idx[1], idx[2]])
        }
        
        return(out)
      }
      outputStat <- mapply(FUN = outOrg, pwMatListOut, statMean,
                           bcStatMean, lowCI, upCI, bcLowCI, bcHighCI,  
                           MoreArgs = list(pw = pw, pwNms = pw_nms),
                           SIMPLIFY = FALSE)
      
      pw_res <- outputStat
      if(fst){
        names(pw_res) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
      } else {
        names(pw_res) <- c("gstEst", "gstEstHed", "djostEst")
      }
      
      # define pwWrite for output
      sprt <- lapply(names(pw_res), FUN = `c`, c("", "", "", "", "", "", ""))
      pwWrite <- lapply(pw_res, function(x){
        comparison <- rownames(x)
        cols <- colnames(x)
        rownames(x) <- NULL
        out <- cbind(comparison, round(x, 4))
        out <- rbind(colnames(out), out)
        colnames(out) <- NULL
        return(out)
      })
      pwWrite <- mapply(FUN = "rbind", sprt, pwWrite, SIMPLIFY = FALSE)
      pwWrite <- do.call("rbind", pwWrite)
      # write results
      if(!is.null(on)){
        if(write_res==TRUE){
          xlsx::write.xlsx(pwWrite, file = paste(of, "[fastDivPart].xlsx", 
                                                 sep = ""),
                           sheetName = "Pairwise_bootstrap", 
                           col.names = FALSE,
                           row.names = FALSE, append = TRUE)
        } else {
          # text file alternatives
          pw_bts <- file(paste(of, "Pairwise-bootstrap[fastDivPart].txt", 
                               sep = ""), "w")
          #cat(paste(colnames(pw_bs_out),sep=""),"\n",sep="\t",file=pw_bts)
          for(i in 1:nrow(pwWrite)){
            cat(pwWrite[i,], "\n", file = pw_bts, sep = "\t")
          }
          close(pw_bts)
        }
      } 
    }
    zzz<-gc()
    rm(zzz)
    ############################################################################
    #pw plotter
    if(plot_res==TRUE && plt==TRUE && bspw==TRUE){
      pwso <- list()
      for(i in 1:length(pw_res)){
        pwso[[i]] <- order(pw_res[[i]][, 1], decreasing = FALSE)
        #if(length(pwso[[i]]) >= 100){
        #  pwso[[i]]<-pwso[[i]][(length(pwso[[i]])-99):length(pwso[[i]])]
        #}
      }
      if(fst){
        names(pwso) <- namer[-c(1:3, length(namer))]
      } else {
        names(pwso) <- namer[-(1:3)]
      }
      
      # define plot parameters 
      plot.call_pw <- list()
      plot.extras_pw <- list()
      xy.labels_pw <- list()
      y.pos_pw <- list()
      x.pos_pw = 1:length(pwso[[i]])
      fn_pre_pw <- list()
      direct = of
      #Plot Gst_Nei
      plot.call_pw[[1]]=c("plot(pw_res[[1]][pwso[[1]],1],
                          ylim=c(0,(max(pw_res[[1]][,3])+
                          min(pw_res[[1]][,3]))),xaxt='n',
                          ylab=names(pw_res)[1],type='n',
                          xlab='Pairwise comparisons 
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[1]]=c("points(pw_res[[1]][pwso[[1]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[1]]),pw_res[[1]][pwso[[1]],6],
                            1:length(pwso[[1]]),pw_res[[1]][pwso[[1]],6],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[5]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[1]] = data.frame(pairwise_name = pw_nms[pwso[[1]]],
                                     Gst_Nei = round(pw_res[[1]][pwso[[1]], 1],4),
                                     Gst_Hedrick = round(pw_res[[2]][pwso[[1]], 1],4),
                                     D_jost = round(pw_res[[3]][pwso[[1]], 1],4))
      
      y.pos_pw[[1]] = pw_res[[1]][pwso[[1]], 1]
      fn_pre_pw[[1]] <- names(pw_res)[1]
      
      
      
      # Plot Gst_Hedrick
      plot.call_pw[[2]]=c("plot(pw_res[[2]][pwso[[2]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(pw_res)[2],type='n',
                          xlab='Pairwise comparisons
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[2]]=c("points(pw_res[[2]][pwso[[2]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[2]]),pw_res[[2]][pwso[[2]],6],
                            1:length(pwso[[2]]),pw_res[[2]][pwso[[2]],7],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[6]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[2]] = data.frame(pairwise_name = pw_nms[pwso[[2]]],
                                     Gst_Nei = round(pw_res[[1]][pwso[[2]],1],4),
                                     Gst_Hedrick = round(pw_res[[2]][pwso[[2]],1],4),
                                     D_jost = round(pw_res[[3]][pwso[[2]],1],4))
      
      y.pos_pw[[2]] = pw_res[[2]][pwso[[2]],1]
      fn_pre_pw[[2]] <- names(pw_res)[2]
      
      
      # Plot D_jost
      plot.call_pw[[3]]=c("plot(pw_res[[3]][pwso[[3]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(pw_res)[3],type='n',
                          xlab='Pairwise comparisons 
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[3]]=c("points(pw_res[[3]][pwso[[3]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[3]]),pw_res[[3]][pwso[[3]],6],
                            1:length(pwso[[3]]),pw_res[[3]][pwso[[3]],7],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[7]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[3]]=data.frame(pairwise_name=pw_nms[pwso[[3]]],
                                   Gst_Nei=round(pw_res[[1]][pwso[[3]],1],4),
                                   Gst_Hedrick=round(pw_res[[2]][pwso[[3]],1],4),
                                   D_jost=round(pw_res[[3]][pwso[[3]],1],4))
      
      y.pos_pw[[3]]=pw_res[[3]][pwso[[3]],1]
      fn_pre_pw[[3]]<-names(pw_res)[3]
      #plot(Fst_WC)
      if(fst==TRUE){
        plot.call_pw[[4]]=c("plot(pw_res[[4]][pwso[[4]],1],
                            ylim=c(0,(max(pw_res[[4]][,3])+
                            min(pw_res[[4]][,3]))),xaxt='n',ylab=names(pw_res)[4],type='n',
                            xlab='Pairwise comparisons 
                            \n (Hover over a point to see pairwise info.)',
                            cex.lab=1.2,cex.axis=1.3,las=1)")
        
        plot.extras_pw[[4]]=c("points(pw_res[[4]][pwso[[4]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],6],
                              1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],7],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=as.numeric(plot_data321[7]),
                              lwd=1,lty=2,col='red')")
        
        xy.labels_pw[[4]]=data.frame(pairwise_name=pw_nms[pwso[[4]]],
                                     Gst_Nei=round(pw_res[[1]][pwso[[4]],1],4),
                                     Gst_Hedrick=round(pw_res[[2]][pwso[[4]],1],4),
                                     D_jost=round(pw_res[[3]][pwso[[4]],1],4),
                                     Fst_WC=round(pw_res[[4]][pwso[[4]],1],4))
        
        y.pos_pw[[4]]=pw_res[[4]][pwso[[4]],1]
        fn_pre_pw[[4]]<-names(pw_res)[4]
      }
    }
    ############################### Bootstrap end ##############################
    
    
    ################################# Plot resuts ###############################
    #make necessary data available
    if(plt==TRUE && plot_res==TRUE && bsls==TRUE && bspw==TRUE){
      pl<-list(bs_res=bs_res,
               pw_res=pw_res,
               accDat=accDat,
               lso123=lso123,
               pwso=pwso,
               plot.call_loci=plot.call_loci,
               plot.extras_loci=plot.extras_loci,
               xy.labels_loci=xy.labels_loci,
               x.pos_loci=x.pos_loci,
               y.pos_loci=y.pos_loci,
               fn_pre_loci=fn_pre_loci,
               direct=direct,
               plot_loci="TRUE",
               plot_pw="TRUE",
               plot.call_pw=plot.call_pw,
               plot.extras_pw=plot.extras_pw,
               xy.labels_pw=xy.labels_pw,
               y.pos_pw=y.pos_pw,
               fn_pre_pw=fn_pre_pw,
               x.pos_pw=x.pos_pw,
               pw=pw,
               plot_data321=plot_data321,
               fst=fst)
    } else if (plt==TRUE && plot_res==TRUE && bsls==TRUE && bspw==FALSE){
      pl<-list(bs_res=bs_res,
               accDat=accDat,
               lso123=lso123,
               plot.call_loci=plot.call_loci,
               plot.extras_loci=plot.extras_loci,
               xy.labels_loci=xy.labels_loci,
               x.pos_loci=x.pos_loci,
               y.pos_loci=y.pos_loci,
               fn_pre_loci=fn_pre_loci,
               direct=direct,
               plot_loci="TRUE",
               plot_pw="FALSE",
               plot_data321=plot_data321,
               fst=fst)
    } else if (plt==TRUE && plot_res==TRUE && bsls==FALSE && bspw==TRUE){
      pl<-list(pw_res=pw_res,
               accDat=accDat,
               pwso=pwso,
               plot.call_pw=plot.call_pw,
               plot.extras_pw=plot.extras_pw,
               xy.labels_pw=xy.labels_pw,
               x.pos_pw=x.pos_pw,
               y.pos_pw=y.pos_pw,
               fn_pre_pw=fn_pre_pw,
               direct=direct,
               plot_loci="FALSE",
               plot_pw="TRUE",
               pw=pw,plot_data321=plot_data321,
               fst=fst)
    }
    if(!is.null(on)){
      if (plt==TRUE && plot_res==TRUE){
        suppressWarnings(plotter(x=pl,img="1000x600"))
      }
    }
    zzz<-gc()
    rm(zzz)
    
    if(pWise | bspw){
      # Create mean pairwise values (for Erin Landguth 12/12)
      meanPairwise <- lapply(pwMatListOut, function(x){
        mean(x, na.rm = TRUE)
      })
      names(meanPairwise) <- names(pwMatListOut)
    }
    
    
    ###########################################################################
    #Data for output
    if(bspw == TRUE && bsls == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_locus = bs_res1,
           bs_pairwise = pw_res,
           bs_pairwise_loci = pwLocOutput)
    } else if(bspw == TRUE && bsls == FALSE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_pairwise = pw_res,
           bs_pairwise_loci = pwLocOutput)
    } else if(bspw == FALSE && bsls == TRUE && pWise == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_locus = bs_res1,
           pw_locus = locstats)
    } else if(bspw == FALSE && bsls == FALSE && pWise == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           pw_locus = locstats)
    } else if(bspw == FALSE && bsls == TRUE && pWise == FALSE){
      list(standard = ot1out,
           estimate = ot2out,
           bs_locus = bs_res1)
    } else if(bspw == FALSE && bsls == FALSE && pWise == FALSE){
      list(standard = ot1out,
           estimate = ot2out)
    }
  }
}
################################################################################
# fastsDivPart end                                                             #
################################################################################