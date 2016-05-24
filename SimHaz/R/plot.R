###--------------------###
### Plotting function  ###
###--------------------###

# functions for incidence plot
#' @export
plot_simuData <- function(data, title="Sample Survival Data"){
  plot.new()
  plot.window(xlim=c(0,max(data$stop)), ylim=c(1,length(unique(data$id))))
  axis(1)
  axis(2)
  #  plot(x=NULL,y=NULL, xlim=c(0, 24), ylim=c(1, nrow(data)))
  title(main=title)
  title(xlab="time")
  title(ylab="Subjects")
  df_new <- data.frame(matrix(ncol = 4, nrow = length(unique(data$id))))
  colnames(df_new) = c('id', "futime",'x',"order")
  df_new$id = unique(data$id)
  id_exp = data[data$x==1, 'id']
  df_new$x = ifelse(df_new$id %in% id_exp, 1, 0)
  for (i in unique(data$id)){
    df_new[df_new$id==i, "futime"] = sum(data[data$id==i, "stop"] - data[data$id==i, "start"])
  }
  df_new <- df_new[order(-df_new$x, -df_new$futime),]
  df_new$order <- seq(1, length(unique(data$id)))
  for (i in 1:nrow(data)) {
    data[i,"index"] = df_new[df_new$id==data[i,"id"],'order']
  }
  points(data[data$status==1, 'stop'], data[data$status==1, 'index'], pch=4)
  points(data[data$status==0, 'stop'], data[data$status==0, 'index'])
  colors = ifelse(data$x==1,'red','blue')
  segments(data$start, data$index, data$stop, data$index, col = colors )
}

plot_simuData_jasa <- function(data){
  plot.new()
  plot.window(xlim=c(0,max(data$tstop)), ylim=c(1,length(unique(data$id))))
  axis(1)
  axis(2)
  #  plot(x=NULL,y=NULL, xlim=c(0, 24), ylim=c(1, nrow(data)))
  title(main="Jasa Data")
  title(xlab="time")
  title(ylab="Subjects")
  df_new <- data.frame(matrix(ncol = 4, nrow = length(unique(data$id))))
  colnames(df_new) = c('id', "futime",'x',"order")
  df_new$id = unique(data$id)
  id_exp = data[data$trans==1, 'id']
  df_new$trans = ifelse(df_new$id %in% id_exp, 1, 0)
  for (i in unique(data$id)){
    df_new[df_new$id==i, "futime"] = sum(data[data$id==i, "tstop"] - data[data$id==i, "tstart"])
  }
  df_new <- df_new[order(-df_new$trans, -df_new$futime),]
  df_new$order <- seq(1, length(unique(data$id)))
  for (i in 1:nrow(data)) {
    data[i,"index"] = df_new[df_new$id==data[i,"id"],'order']
  }
  colors = ifelse(data$trans==1,'red','blue')
  points(data[data$death==1, 'tstop'], data[data$death==1, 'index'], pch=4)
  points(data[data$death==0, 'tstop'], data[data$death==0, 'index'])
  segments(data$tstart, data$index, data$tstop, data$index, col = colors )
}

# function to plot power curves
#' @export

plot_power<-function(table_df,N,type,exp.prop,min.futime,min.postexp.futime,show.plot=FALSE,newplot=FALSE,
                     col=NULL,lty=NULL,lwd=NULL,pch=NULL){

  df<-subset(table_df, table_df$i_N==N & table_df$i_type==type & table_df$i_exp.prop == exp.prop & 
               table_df$i_min.futime==min.futime & table_df$i_min.postexp.futime==min.postexp.futime,
             select=c('i_N','N_eff','i_beta','pow'))
  if(show.plot){
    if(newplot){
      plot(unique(df$i_beta),df[seq(1,nrow(df),nrow(df)/length(unique(df$i_beta))),]$pow,main="Random censoring time",xlab=expression(beta),ylab="Statistical power",
           type="b",col=col,lty=lty,lwd=lwd,ylim=c(0,1),pch=pch)
    }
    else{
      lines(unique(df$i_beta),df[seq(1,nrow(df),nrow(df)/length(unique(df$i_beta))),]$pow,type="b",col=col,lty=lty,lwd=lwd,pch=pch)
    }
  }
  return(df)
}


