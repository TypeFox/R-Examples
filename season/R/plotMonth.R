# plotMonth.R
# function to plot monthly results by month
# assumes the data are in order
# assumes there are numeric variables for month and year

plotMonth<-function(data,resp,panels=12, ...){

  year <- yaxis <- Month <- NULL # Setting some variables to NULL first (for R CMD check)
  
  if (panels!=1&panels!=12){stop("panels must be 1 or 12")}
  data$yaxis=subset(data,select=resp)[,1]
  
  # 12 panels
  data$Month=factor(data$month,levels=1:12,labels=month.abb) # to change facet_wrap labels
  if(panels==12){
  gplot = ggplot(data, aes(year, yaxis)) +
    geom_line()+
    theme_bw()+
    xlab(' ') +        
    ylab(resp) +        
    facet_wrap(~Month)
  print(gplot)
  }
  
  # 1 panels
  if(panels==1){
    gplot = ggplot(data, aes(year, yaxis,color=Month)) +
      geom_line()+
      theme_bw()+
      xlab(' ') +        
      ylab(resp)         
    print(gplot)
  }
  
}
