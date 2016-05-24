getColorsForGroups <- function(group,colors=c("red","blue","green","purple","grey","cyan","brown","pink")){
  cluster_colors=group
  
  if(max(group)<=length(colors)){
    for(i in 1:max(group)){
      cluster_colors[which(group==i)] <- colors[i]
    }
    return(cluster_colors)
  }
  else{
    print("ERROR: Not enough colors using the default color argument for the different groups, PLEASE inform the colors argument")
    return(NULL)
  }
}