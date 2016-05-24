plotORdensity <-
function(x, log = "x", ncol = 2, adjust.zeroinf = TRUE, zero.pos = 0.005,
                           inf.pos = 200, binwidth = 0.15, color = "black", xlab = "odds ratio",
                           save = FALSE, path = NULL, file.name = NULL,
                           format = "jpeg",...){
  
  
  ##### internal functions ############
  plot.ORdensity.ampliconduo <-
    function(x, log, ncol, adjust.zeroinf, zero.pos,
             inf.pos, binwidth, color, xlab,
             save, path, file.name, format, main = NULL, ...){
      if(is.null(main)){
        main = x[1,9]
      }
      
      OR <- ..density.. <- NULL ## just for CMD check not to complain
      
      if(adjust.zeroinf == T){
        d2 <- x
        d2[d2[4] == 0,][4] <- zero.pos ### d2[4] is same as d2$OR
        d2[d2[4] == Inf,][4] <- inf.pos
        plotOR <-qplot(
          OR,
          log = log,
          data = d2,
          ..density..,
          geom = "histogram",
          main = main,
          fill = I(color),
          binwidth = binwidth,
          ...
        ) +
          xlab(xlab)
        
      }else{
        plotOR <- qplot(
          OR,
          log = log,
          data = x,
          ..density..,
          geom = "histogram",
          main = main,
          fill = I(color),
          binwidth = binwidth,
          ...
        ) +
          xlab(xlab)
      } 
      
      if(save == T){
        if(is.null(file.name)){
          file.name = paste("ORdensity", "_", Sys.Date(), ".",format , sep ="")
        }else{
          file.name = paste(file.name, ".", format, sep = "")
        }
        ggsave(filename = file.name, path = path, ... )
      }
      print(plotOR) 
    }
  
  plot.ORdensity.ampliconduo.set <-
    function(x, log, ncol, adjust.zeroinf, zero.pos,
             inf.pos, binwidth, color, xlab,
             save, path, file.name, format, ...){
      
      OR <- NULL ## just to make CMD check not complain
      ..density.. <- NULL
      
      if(length(x)> 1){
        dat = x[[1]]
        for(i in 2:length(x)){
          dat = rbind(dat, x[[i]])
        }
      }else{
        dat = x[[1]]
      }
      
      if(adjust.zeroinf == T){
        d2 <- dat
        d2[d2[4] == 0,][4] <- zero.pos ## d2[4] is same as d2$OR
        d2[d2[4] == Inf,][4] <- inf.pos
        plotOR <- qplot(
          OR,
          log = log,
          data = d2,
          ..density..,
          geom = "histogram",
          binwidth = binwidth,
          fill = I(color),
          ...
        ) +
          xlab(xlab) +
          facet_wrap(~sample, ncol = ncol)
        
        
      }else{
        plotOR <- qplot(
          OR,
          log = log,
          data = dat,
          ..density..,
          geom = "histogram",
          binwidth = binwidth,
          fill = I(color),
          ...
        ) +
          xlab(xlab) +
          facet_wrap(~sample, ncol = ncol)
      }
      
      if(save == T){
        if(is.null(file.name)){
          file.name = paste("ORdensity", "_", Sys.Date(), ".", format, sep ="")
        }else{
          file.name = paste(file.name, ".", format, sep = "")
        }
        ggsave(filename = file.name, path = path, ... )
      }
      print(plotOR)
    }
  
  ##################################
  ncol=as.integer(ncol)
  xclass <- class(x)
  if(xclass[1] == "list"){
    plot.ORdensity.ampliconduo.set(x, log, ncol, adjust.zeroinf, zero.pos,
                                   inf.pos, binwidth, color, xlab,
                                   save, path, file.name, format, ...)
  }else{
    if(xclass[1] == "data.frame"){
      plot.ORdensity.ampliconduo(x, log, ncol, adjust.zeroinf, zero.pos,
                                 inf.pos, binwidth, color, xlab,
                                 save, path, file.name, format, ...)
    }else{
      stop("wrong data format of x!")
    }
  } 
}
