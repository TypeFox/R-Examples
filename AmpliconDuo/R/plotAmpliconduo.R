plotAmpliconduo <-
function(x , color.treshold = 0.05, 
                             xlab = "Abundance (PCR A)",
                             ylab = "Abundance (PCR B)",
                             main = NULL, log = "xy",
                             corrected = TRUE, asp = 1, 
                             legend.position = NULL, 
                             save = FALSE, path = NULL, 
                             file.name = NULL, format = "jpeg",
                             h.start = 0, ...){
  
  if(is.null(main)){
    main = x[1,9]
  }
  
  freqA <- freqB <-  q <- p <- NULL ## just for CMD check not to complain

  if(corrected == T){
    legend.title = paste("q > ", color.treshold, sep = "")
    ampliplot <- qplot(
      freqA,
      freqB,
      xlab = xlab,
      ylab = xlab,
      data=x,
      log=log,
      asp = asp,
      colour = q > color.treshold,
      main = main, ...
    ) + scale_colour_discrete(name = legend.title,  h.start = h.start)
  }else{
    legend.title = paste("p > ", color.treshold, sep = "")
    ampliplot <- qplot(
      freqA,
      freqB,
      xlab = xlab,
      ylab = xlab,
      data=x,
      log=log,
      asp = asp,
      colour = p > color.treshold,
      main = main, ...
    ) + scale_colour_discrete(name = legend.title,  h.start = h.start)
    
  }
  if(save == T){
    if(is.null(file.name)){
      file.name = paste(main, "_", Sys.Date(), ".", format , sep ="")
    }else{
      file.name = paste(file.name, ".", format, sep = "")
    }
    ggsave(filename = file.name, path = path, ... )
  }  
  print(ampliplot)
}
