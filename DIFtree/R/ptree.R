ptree <-
function(info,
                  item,
                  params,
                  X,
                  cex.lines,
                  cex.main,
                  cex.branches,
                  cex.coefs,
                  title){
  
  # compute number of nodes 
  endnodes <- get_endnodes(info)
  
  n_levels <- length(unique(info[,"level"]))
  
  hilfspunkte <- list()
  hilfspunkte[[1]] <- matrix(NA,nrow=2^n_levels,ncol=2)
  hilfspunkte[[1]][,1] <- 2^n_levels
  hilfspunkte[[1]][,2] <- rep(n_levels+1,2^n_levels)
  
  steps <- 2^((n_levels:1-1))
  
  for(i in 1:n_levels){
    
    hilfspunkte[[i+1]] <- hilfspunkte[[i]]
    hilfspunkte[[i+1]][,2] <- rep(n_levels+1-i,2^n_levels)
    
    help  <- c(-steps[i],steps[i])
    help1 <- rep(help,each=steps[i])
    help2 <- rep(help1,length=2^n_levels) 
    hilfspunkte[[i+1]][,1] <- hilfspunkte[[i]][,1]+help2 
    
    which_knots <- info[info[,"level"]==i,"node"]
    help3 <- seq(1,2^n_levels)
    help4 <- split(help3,rep(1:2^(i-1),each=2^n_levels/2^(i-1)))
    help5 <- unlist(lapply(which_knots, function(j) help4[[j]]))
    hilfspunkte[[i+1]][-help5,] <- hilfspunkte[[i]][-help5,]
    
  }
  
  
  plot.new()
  plot.window(ylim=c(0.5,n_levels+1),xlim=c(0,2^(n_levels+1)))
  rect(0,0.5,2^(n_levels+1),n_levels+1, border = grey(0.9),col = grey(0.9))
  
  
  for(j in length(hilfspunkte):2){
    for(i in 1:(2^n_levels)){
      lines(c(hilfspunkte[[j-1]][i,1],hilfspunkte[[j]][i,1]),c(hilfspunkte[[j-1]][i,2],hilfspunkte[[j]][i,2]),
            lwd=cex.lines)
    }
  }
  if(is.null(title)){
    title <- paste("Item",item)
  }
  title(title,cex.main=cex.main)
  
  # add estimates
  params_hat <- round(params,3)
  points_params <- unique(hilfspunkte[[n_levels+1]])
  points_params[,2] <- points_params[,2]-0.2
  for(i in 1:length(params_hat)){
    draw.ellipse(x=points_params[i,1],y=points_params[i,2],a=0.8,b=0.2,lwd=cex.lines,col=grey(0.8))
  }
  text(points_params[,1],points_params[,2],params_hat,cex=cex.coefs)
  
  # add numbers
  points_params[,2] <- points_params[,2]+0.2
  for(i in 1:length(params_hat)){
    rect(points_params[i,1]-0.2,points_params[i,2]-0.1,points_params[i,1]+0.2,points_params[i,2]+0.1,col=grey(0.9),lwd=cex.lines)
    text(points_params[i,1],points_params[i,2],endnodes[i],cex=cex.branches)
  }
  
  # add infos 
  for(i in 1:nrow(info)){
    help4 <- split(help3,rep(1:2^(info[i,"level"]-1),each=2^n_levels/2^(info[i,"level"]-1)))[[info[i,"node"]]]
    point_var <- unique(hilfspunkte[[info[i,"level"]]][help4,])
    points(point_var[1],point_var[2],cex=cex.lines-1,pch=19)
    point_left  <- c(point_var[1]-steps[info[i,"level"]],point_var[2]-0.5)
    point_right <- c(point_var[1]+steps[info[i,"level"]],point_var[2]-0.5)
    var   <- info[i,"variable"]
    thres <- info[i,"threshold"]
    sort_values <- unique(sort(X[,var]))
    if(thres==min(sort_values)){
      text(point_left[1],point_left[2],paste0(var,"=",round(thres,2)),cex=cex.branches,adj=c(1,0))
    } else{
      text(point_left[1],point_left[2],paste0(var,"<=",round(thres,2)),cex=cex.branches,adj=c(1,0))
    }
    if(thres==max(sort_values[-length(sort_values)])){
      text(point_right[1],point_right[2],paste0(var,"=",round(max(sort_values),2)),cex=cex.branches,adj=c(0,0))
    } else{
      text(point_right[1],point_right[2],paste0(var,">",round(thres,2)),cex=cex.branches,adj=c(0,0))
    }
  }
  
}
