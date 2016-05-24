plotAmpliconduo.set <-
function(x , color.treshold = 0.05,
                                 xlab = "Abundance (PCR A)",
                                 ylab = "Abundance (PCR B)", 
                                 log = "xy",  corrected = TRUE, 
                                 asp = 1, nrow = 1, legend.position = NULL ,
                                 save = FALSE, path = NULL, file.name = NULL, 
                                 format = "jpeg", h.start = 0, ...){
  
  if(length(x)> 1){
    dat = x[[1]]
    for(i in 2:length(x)){
      dat = rbind(dat,x[[i]])
    }
  }else{
    dat =x[[1]]
  }
  
  freqA <- freqB <- sample <- q <- p <- NULL ## just for CMD check not to complain
  
  ## determines legend position
  if(is.null(legend.position)){
    s.z <- length(x)
    rest <- s.z %% nrow 
    mod <- s.z %/% nrow
    if (rest == 0 || rest > s.z || (mod == 1 && rest > 0)){
      if(mod > nrow){
        legend.position <- "top"
      }else{
        legend.position <- "right"
      }
    }else{
      legend.position <- c(1, 0)
    }
  }
  
  if(corrected == TRUE){
    legend.title = paste("q > ", color.treshold, sep = "")
    ampliplot <- qplot(
      freqA,
      freqB,
      xlab = xlab,
      ylab = ylab,
      data = dat,
      log = log,
      colour= q > color.treshold,
      asp = asp,
      ...
    ) + 
      facet_wrap(~ sample, nrow = nrow) +
      theme(legend.position = legend.position, legend.justification = legend.position) + 
      scale_colour_discrete(name = legend.title, h.start = h.start)
    
  }else{
    legend.title = paste("p > ", color.treshold, sep = "")
    ampliplot <- qplot(
      freqA,
      freqB,
      xlab = xlab,
      ylab = ylab,
      data = dat,
      log = log,
      colour= p > color.treshold,
      asp = asp,
      ...
    ) + 
      facet_wrap(facets = ~ sample, nrow = nrow) +
      theme(legend.position = legend.position) + scale_colour_discrete(name = legend.title)
  } 
  
  if(save == T){
    if(is.null(file.name)){
      file.name = paste("ampliconduo_", Sys.Date(), ".", format , sep ="")
    }else{
      file.name = paste(file.name, ".", format, sep = "")
    }
    ggsave(filename = file.name, path= path, ... )
  }
  
  print(ampliplot)
}
