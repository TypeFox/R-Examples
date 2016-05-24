overview <-
function(r,row_col=c(3,3)){
  
  ## dynRB_VPa
  
  if(names(r)[[1]]=="plot_data_niche_size_V1"){
    par(mfrow = row_col,
        oma = c(0, 0, 2, 0),
        mar = c(5, 5, 5, 5))
    k<-1
    h<-1
    
    for(i in 1:sqrt(length(r$plot_data_niche_size_V1))){
      plot(r$plot_alpha_grid[[i]][-length(r$plot_alpha_grid[[i]])],r$plot_data_niche_size_V1[[i]][2:length(r$plot_data_niche_size_V1[[i]])],col=4,type="l",ylim=c(0,1.1),xlab="alpha","ylab"="volume prod")
      title(main = paste("S1: ",names(r$plot_data_niche_size_V1)[[i]],"\nNiche Size: Integral=", round(r$plot_data_niche_size_V1[[i]][1],4)))
      if(k/h==row_col[1]*row_col[2]){
        title("Niche size of species S1", outer = T)
        h<-h+1
        readline(prompt = "Press <Enter> to continue...")}
      k<-k+1
    }
    title("Niche size of species S1", outer = T)
    
    readline(prompt = "Press <Enter> to continue...")
    
    par(mfrow = row_col,
        oma = c(0, 0, 2, 0),
        mar = c(5, 5, 5, 5))
    k<-1
    h<-1
    for(i in 1:length(r$plot_data_overlap)){
      plot(r$plot_alpha_grid[[i]],r$plot_data_overlap[[i]][2:length(r$plot_data_overlap[[i]])],col=4,type="l",ylim=c(0,1.1),xlab="alpha","ylab"="volume prod")
      title(main = paste("S1: ",names(r$plot_data_niche_size_V1)[[i]],"\nS2: ",names(r$plot_data_niche_size_V2)[[i]],"\nOverlap:  Integral=", round(r$plot_data_overlap[[i]][1],4)))
      if(k/h==row_col[1]*row_col[2]){
        title("overlap S1 in S2 (product)", outer = T)
        h<-h+1}
      k<-k+1
    }
    title("overlap S1 in S2 (product)", outer = T)
  }
  
  ## dynRB_Pn 
  
  if(names(r)[[1]]=="names_data_V1"){
    names_data<-names(as.data.frame(r$plot_data_overlap[[1]]))
    par(mfrow = row_col,
        oma = c(0, 0, 2, 0),
        mar = c(5, 5, 5, 5))
    h<-1
    k<-1
    plot_alpha_grid<-seq(0,1,length=length(r$plot_data_overlap[[1]][2:length(r$plot_data_overlap[[1]][,1]),1]))
    for(i in 1:length(r$plot_alpha_grid)){
      for(j in 1:ncol(r$plot_data_overlap[[i]])){
        plot(plot_alpha_grid,r$plot_data_overlap[[i]][2:length(r$plot_data_overlap[[i]][,1]),j],type="l",ylim=c(0,1.1),xlab="alpha","ylab"="overlap",col=4,)             
        title(main = paste("S1: ",r$names_data_V1[[i]],"\nS2: ",r$names_data_V2[[i]],"\n",names_data[j],"Overlap=", round(r$plot_data_overlap[[i]][1,j],4)))
        if(k/h==row_col[1]*row_col[2]){
          title("overlap S1 in S2 (product) per dimension", outer = T)
          h<-h+1
          readline(prompt = "Press <Enter> to continue...")}
        k<-k+1
      }
    }
    title("overlap S1 in S2 (product) per dimension", outer = T)
  }
  
  ## dynRB_Vn 
  
  if(names(r)[[1]]=="plot_data_volume"){
    par(mfrow = row_col,
        oma = c(0, 0, 2, 0),
        mar = c(5, 5, 5, 5))
    h<-1
    k<-1
    for(i in 1:length(r$result[,1])){
      name<-r$result[i,1]
      for(j in 1:(length(r$result[1,])-1)){
        plot(r$plot_alpha_grid[-length(r$plot_alpha_grid)],r$plot_data_volume[[i]][,j],type="l",col=4,ylim=c(0,1.1),xlab="alpha","ylab"="volume")
        title(main = paste("S1: ",name,",\n",r$names_data_S1[j],": Integral=", round(r$result[i,j+1],4)))
        if(k/h==row_col[1]*row_col[2]){
          title("Niche size of species S1 per dimension", outer = T)
          h<-h+1
          readline(prompt = "Press <Enter> to continue...")}
        k<-k+1
      }
    }
    title("Niche size of species S1 per dimension", outer = T)
  }
  
}
