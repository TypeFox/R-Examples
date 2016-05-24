################################################################################
# plotter, a function to create interactive plots of results from divPart      #
################################################################################
plotter<-function(x,img="1200x600"){
  if ("sendplot" %in% rownames(installed.packages())) {
    requireNamespace("sendplot")
    x=x
    spot.radius=5
    jjj<-x
    fl_ext<-c(".tif","Dot.png","Dot.tif")
    bs_res<-list()
    lso123<-list()
    accDat<-list()
    sp.header<-list()
    pw_res<-list()
    pwso<-list()
    pw<-list()
    plot_data321<-list()
    if(jjj$plot_loci==TRUE && jjj$plot_pw==FALSE){
      bs_res<<-jjj$bs_res
      lso123<<-jjj$lso123
      accDat<<-jjj$accDat
      if(length(lso123) > 150){
        image.size <- "2400x1200"
      } else {
        image.size=img
      }
      #Gst_loci
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[1]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[1]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[1]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[1]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[1]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      #clean up
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[1]],"_locus_stat_",fl_ext,sep=""))
      #G'st_loci
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[2]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[2]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[2]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[2]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[2]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[2]],"_locus_stat_",fl_ext,sep=""))
      #Djost_loci
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[3]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[3]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[3]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[3]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[3]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[3]],"_locus_stat_",fl_ext,sep=""))
      #Fst_WC_loci
      if(jjj$fst==TRUE){
        suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[4]],
                                   x.pos=jjj$x.pos_loci,
                                   y.pos=jjj$y.pos_loci[[4]],
                                   xy.type="points",
                                   plot.extras=jjj$plot.extras_loci[[4]],
                                   mai.mat=NA,
                                   mai.prc=FALSE,
                                   xy.labels=jjj$xy.labels_loci[[4]],
                                   image.size=image.size,
                                   spot.radius=5,
                                   fname.root=paste(jjj$fn_pre_loci[[4]],
                                                    "_locus_stat_",sep=""),
                                   dir=jjj$direct,
                                   window.size="2100x1000"))
        unlink(paste(jjj$direct,jjj$fn_pre_loci[[4]],"_locus_stat_",fl_ext,sep=""))
      }
      if(exists("jjj", where=".GlobalEnv")==TRUE){
        rm(jjj, pos=".GlobalEnv")
      }
      if(exists("accDat", where=".GlobalEnv")==TRUE){
        rm(accDat, pos=".GlobalEnv")
      }
      if(exists("bs_res", where=".GlobalEnv")==TRUE){
        rm(bs_res, pos=".GlobalEnv")
      }
      if(exists("lso123", where=".GlobalEnv")==TRUE){
        rm(lso123, pos=".GlobalEnv")
      }
      if(exists("sp.header", where=".GlobalEnv")==TRUE){
        rm(sp.header, pos=".GlobalEnv")
      }
      if(exists("plot_data321", where=".GlobalEnv")==TRUE){
        rm(plot_data321, pos=".GlobalEnv")
      }
      #rm(jjj,accDat,bs_res,lso123,sp.header,pos=".GlobalEnv")
      
    } else if(jjj$plot_loci==FALSE && jjj$plot_pw==TRUE){
      accDat<<-jjj$accDat
      pw_res<<-jjj$pw_res
      pwso<<-jjj$pwso
      pw<<-jjj$pw
      plot_data321<<-jjj$plot_data321
      if(length(pwso) > 150){
        image.size <- "2400x1200"
      } else {
        image.size=img
      }
      #Gst_pw
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[1]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[1]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[1]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[1]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[1]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[1]],"_pairwise_stats_",fl_ext,sep=""))
      #G'st_pw
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[2]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[2]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[2]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[2]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[2]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[2]],"_pairwise_stats_",fl_ext,sep=""))
      #Djost_pw
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[3]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[3]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[3]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[3]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[3]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",fl_ext,sep=""))
      #Fst_WC_pw
      if(jjj$fst==TRUE){
        suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[4]],
                                   x.pos=jjj$x.pos_pw,
                                   y.pos=jjj$y.pos_pw[[4]],
                                   xy.type="points",
                                   plot.extras=jjj$plot.extras_pw[[4]],
                                   mai.mat=NA,
                                   mai.prc=FALSE,
                                   xy.labels=jjj$xy.labels_pw[[4]],
                                   image.size=image.size,
                                   spot.radius=5,
                                   fname.root=paste(jjj$fn_pre_pw[[4]],
                                                    "_pairwise_stats_",sep=""),
                                   dir=jjj$direct,
                                   window.size="2100x1000"))
        unlink(paste(jjj$direct,jjj$fn_pre_pw[[4]],"_pairwise_stats_",
                     fl_ext,sep=""))
      }
      
      if(exists("jjj", where=".GlobalEnv")==TRUE){
        rm(jjj, pos=".GlobalEnv")
      }
      if(exists("accDat", where=".GlobalEnv")==TRUE){
        rm(accDat, pos=".GlobalEnv")
      }
      if(exists("pw_res", where=".GlobalEnv")==TRUE){
        rm(pw_res, pos=".GlobalEnv")
      }
      if(exists("pwso", where=".GlobalEnv")==TRUE){
        rm(pwso, pos=".GlobalEnv")
      }
      if(exists("sp.header", where=".GlobalEnv")==TRUE){
        rm(sp.header, pos=".GlobalEnv")
      }
      if(exists("plot_data321", where=".GlobalEnv")==TRUE){
        rm(plot_data321, pos=".GlobalEnv")
      }
      if(exists("pw", where=".GlobalEnv")==TRUE){
        rm(pw, pos=".GlobalEnv")
      }
      #rm(jjj,accDat,plot_data,pw,pw_res,pwso,sp.header,pos=".GlobalEnv")
      
    } else if(jjj$plot_loci==TRUE && jjj$plot_pw==TRUE){
      bs_res<<-jjj$bs_res
      lso123<<-jjj$lso123
      accDat<<-jjj$accDat
      pw_res<<-jjj$pw_res
      pwso<<-jjj$pwso
      pw<<-jjj$pw
      plot_data321<<-jjj$plot_data321
      if(length(lso123) > 150 && length(pwso) > 150){
        pwimage.size <- "2400x1200"
        locimage.size <- "2400x1200"
      } else if(length(lso123) > 150 && length(pwso) <= 150){
        pwimage.size <- img
        locimage.size <- "2400x1200"
      } else if(length(lso123) <= 150 && length(pwso) > 150){
        pwimage.size <- "2400x1200"
        locimage.size <- img
      } else {
        locimage.size <- img
        pwimage.size <- img
      }
      #Gst_loci
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[1]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[1]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[1]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[1]],
                                 image.size=locimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[1]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[1]],"_locus_stat_",fl_ext,sep=""))
      #G'st_loci
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[2]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[2]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[2]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[2]],
                                 image.size=locimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[2]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[2]],"_locus_stat_",fl_ext,sep=""))
      #Djost_loci
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[3]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[3]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[3]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[3]],
                                 image.size=locimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[3]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[3]],"_locus_stat_",fl_ext,sep=""))
      #Fst_WC_loci
      if(jjj$fst==TRUE){
        suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_loci[[4]],
                                   x.pos=jjj$x.pos_loci,
                                   y.pos=jjj$y.pos_loci[[4]],
                                   xy.type="points",
                                   plot.extras=jjj$plot.extras_loci[[4]],
                                   mai.mat=NA,
                                   mai.prc=FALSE,
                                   xy.labels=jjj$xy.labels_loci[[4]],
                                   image.size=locimage.size,
                                   spot.radius=5,
                                   fname.root=paste(jjj$fn_pre_loci[[4]],
                                                    "_locus_stat_",sep=""),
                                   dir=jjj$direct,
                                   window.size="2100x1000"))
        unlink(paste(jjj$direct,jjj$fn_pre_loci[[4]],"_locus_stat_",fl_ext,sep=""))
      }
      #Gst_pw
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[1]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[1]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[1]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[1]],
                                 image.size=pwimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[1]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[1]],"_pairwise_stats_",
                   fl_ext,sep=""))
      #G'st_pw
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[2]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[2]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[2]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[2]],
                                 image.size=pwimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[2]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[2]],"_pairwise_stats_",fl_ext,sep=""))
      #Djost_pw
      suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[3]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[3]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[3]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[3]],
                                 image.size=pwimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[3]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",
                   fl_ext,sep=""))
      #Fst_WC_pw
      if(jjj$fst==TRUE){
        suppressWarnings(sendplot::imagesend(plot.call=jjj$plot.call_pw[[4]],
                                   x.pos=jjj$x.pos_pw,
                                   y.pos=jjj$y.pos_pw[[4]],
                                   xy.type="points",
                                   plot.extras=jjj$plot.extras_pw[[4]],
                                   mai.mat=NA,
                                   mai.prc=FALSE,
                                   xy.labels=jjj$xy.labels_pw[[4]],
                                   image.size=pwimage.size,
                                   spot.radius=5,
                                   fname.root=paste(jjj$fn_pre_pw[[4]],
                                                    "_pairwise_stats_",sep=""),
                                   dir=jjj$direct,
                                   window.size="2100x1000"))
        unlink(paste(jjj$direct,jjj$fn_pre_pw[[4]],"_pairwise_stats_",
                     fl_ext,sep=""))
      }
      if(exists("jjj", where=".GlobalEnv")==TRUE){
        rm(jjj, pos=".GlobalEnv")
      }
      if(exists("accDat", where=".GlobalEnv")==TRUE){
        rm(accDat, pos=".GlobalEnv")
      }
      if(exists("pw_res", where=".GlobalEnv")==TRUE){
        rm(pw_res, pos=".GlobalEnv")
      }
      if(exists("pwso", where=".GlobalEnv")==TRUE){
        rm(pwso, pos=".GlobalEnv")
      }
      if(exists("sp.header", where=".GlobalEnv")==TRUE){
        rm(sp.header, pos=".GlobalEnv")
      }
      if(exists("plot_data321", where=".GlobalEnv")==TRUE){
        rm(plot_data321, pos=".GlobalEnv")
      }
      if(exists("pw", where=".GlobalEnv")==TRUE){
        rm(pw, pos=".GlobalEnv")
      }
      if(exists("bs_res", where=".GlobalEnv")==TRUE){
        rm(bs_res, pos=".GlobalEnv")
      }
      if(exists("lso123", where=".GlobalEnv")==TRUE){
        rm(lso123, pos=".GlobalEnv")
      }
    } 
  } else {
    stop("The 'sendplot' package must be installed to generate interactive plots.")
  }
}
################################################################################
# plotter end                                                                  #
################################################################################