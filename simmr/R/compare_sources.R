compare_sources = function(simmr_out,source_names=simmr_out$input$source_names,group=1,plot=TRUE) {

# Function to compare between sources within a group both via textual output and with boxplots
# Things to supply are:
# If two sources are given: 
#   - provide the probability of one group having higher dietary proportion than the other
#   - give the probability distribution of the difference
#   - optional boxplot of two 
# If more than two sources are given:
#   - provide the top most likely orderings of the sources
# An optional boxplot of the sources
  
# Throw an error if only one group is specified
if(length(source_names)==1) stop("Use compare_between_groups if you want to compare a single source between groups.")
  
# Throw an error if the source name given doesn't match the source names
if(!all(source_names%in%simmr_out$input$source_names)) stop("Some source names not found in the current source names. Be sure to check case and spelling")
  
# Start with two groups version
if(length(source_names)==2) {
  # Get the output for this particular source on these two groups  
  out_all_src_1 = do.call(rbind,simmr_out$output[[group]])[,source_names[1]]
  out_all_src_2 = do.call(rbind,simmr_out$output[[group]])[,source_names[2]]
  # Produce the difference between the two
  out_diff = out_all_src_1 - out_all_src_2
  cat(paste("Prob ( proportion of",source_names[1],'> proportion of',source_names[2],') =',round(mean(out_diff>0),3)))
  
  if(plot) {
    # Stupid fix for packaging ggplot things
    Source = Proportion = NULL
    df = data.frame(Proportion=c(out_all_src_1,out_all_src_2),Source=c(rep(source_names[1],length(out_all_src_1)),rep(source_names[2],length(out_all_src_2))))
    p = ggplot(df,aes(x=Source,y=Proportion,fill=Source)) + geom_boxplot(alpha=0.5,outlier.size=0) + theme_bw() + theme(legend.position='none') + ggtitle(paste("Comparison of dietary proportions for sources",source_names[1],'and',source_names[2]))
    print(p)
  }
  
} 

# Now for more groups  
if(length(source_names)>2) {
  # Get the output for all the groups
  out_all = do.call(rbind,simmr_out$output[[group]])[,source_names]

  # Now find the ordering of each one
  ordering_num = t(apply(out_all,1,order,decreasing=TRUE))
  Ordering = rep(NA,length=nrow(ordering_num))
  for(i in 1:length(Ordering)) Ordering[i] = paste0(source_names[ordering_num[i,]],collapse=" > ")
  cat('Most popular orderings are as follows:\n')
  tab = t(t(sort(table(Ordering,dnn=NULL),decreasing=TRUE)))
  colnames(tab) = 'Probability'
  # Do not print all of it if too long
  if(nrow(tab)>30) {
    print(round(tab[1:30,]/length(Ordering),4))
  } else {
    print(round(tab/length(Ordering),4))
  }
    
  if(plot) {
    # Stupid fix for packaging ggplot things
    Source = Proportion = NULL
    df = reshape2::melt(out_all)[,2:3]
    colnames(df) = c('Source','Proportion')
    p = ggplot(df,aes(x=Source,y=Proportion,fill=Source)) + 
      scale_fill_viridis(discrete=TRUE) + 
      geom_boxplot(alpha=0.5,outlier.size=0) + 
      theme_bw() + 
      theme(legend.position='none') + 
      ggtitle(paste("Comparison of dietary proportions between sources"))
    print(p)
  }
  
}  

# Return output
if(length(source_names)==2) {
  invisible(list(out_diff))
} else {
  invisible(list(Ordering=Ordering,out_all=out_all))
}  

}