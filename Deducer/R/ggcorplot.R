ggcorplot <- function(cor.mat,data=NULL,lines=TRUE,line.method=c("lm","loess"),type="points",
		alpha=.25,main="auto",var_text_size=5,
		cor_text_limits=c(5,25),level=.05){
	x_var <- y_var <- trans <- rsq <- p <- x_label <- NULL
	#define a helper function (borrowed from the "ez" package)
	ezLev<-function(x,new_order){
		for(i in rev(new_order)){
			x<-relevel(x,ref=i)
		}
		return(x)
	}							
	
	if(all(line.method==c("lm","loess")))
		line.method<-"lm"	
	
	nm <- names(cor.mat)
	for(i in 1:length(nm))
		dat <- if(i==1) d(eval(parse(text=nm[i]),data,parent.frame())) else d(dat, eval(parse(text=nm[i]),data,parent.frame()))
	data <- dat
	names(data) <- nm
	# normalize data
	for(i in 1:length(data)){
		data[,i]<-as.numeric(data[,i])
		data[,i]<-(data[,i]-mean(data[,i],na.rm=TRUE))/sd(data[,i],na.rm=TRUE)
	}
	# obtain new data frame
	z<-data.frame()
	i <- 1
	j <- i
	while(i<=length(data)){
		if(j>length(data)){
			i<-i+1
			j<-i
		}else{
			x <- data[,i]
			y <- data[,j]
			temp<-as.data.frame((cbind(x,y)))
			temp<-cbind(temp,names(data)[i],names(data)[j])
			z<-rbind(z,temp)
			
			j<-j+1
		}
	}
	z<-cbind(z,alpha)
	names(z)=c('x_var','y_var','x_label','y_label','trans')
	z$x_label <- ezLev(factor(z$x_label),names(data))
	z$y_label <- ezLev(factor(z$y_label),names(data))
	z=z[z$x_label!=z$y_label,]
	#obtain correlation values
	z_cor <- data.frame()
	i <- 1
	j <- i
	while(i<=length(data)){
		if(j>length(data)){
			i<-i+1
			j<-i
		}else{
			x <- na.omit(data[,i])
			y <- na.omit(data[,j])
			x_mid <- min(x)+diff(range(x))/2
			y_mid <- min(y)+diff(range(y))/2
			this_cor <- cor.mat[[i]][[j]]$estimate
			this_cor.test <- cor.mat[[i]][[j]]
			this_col <- ifelse(this_cor.test$p.value<level,"red"
					,"blue")
			this_size <- (this_cor)^2
			cor_text <- ifelse(
					this_cor>0
					,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
					,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
			)
			b<-as.data.frame(cor_text)
			b<-cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
			z_cor<-rbind(z_cor,b)
			j<-j+1
		}
	}
	names(z_cor)<-c('cor','x_mid','y_mid','p','rsq','x_label','y_label')
	z_cor$x_label <- ezLev(factor(z_cor$x_label),names(data))
	z_cor$y_label <- ezLev(factor(z_cor$y_label),names(data))
	diag <- z_cor[z_cor$x_label==z_cor$y_label,]
	z_cor<-z_cor[z_cor$x_label!=z_cor$y_label,]
	
	#start creating layers	
	points_layer <- geom_point(aes(x = x_var, y = y_var, alpha=trans),data = z)
	bin_layer<-geom_hex(data = z, mapping = aes(x = x_var, y = y_var,alpha=trans),bins=10)
	lm_line_layer <- stat_smooth(aes(x=x_var, y=y_var), method=line.method)
	cor_text <- geom_text(aes(x=y_mid, y=x_mid, label=cor, size = rsq, colour = p), data=z_cor)
	var_text <- geom_text(aes(x=y_mid, y=x_mid, label=x_label), data=diag, size=var_text_size)

	f <- facet_grid(y_label~x_label,scales='free')
	o <- theme(
			panel.grid.minor = element_blank()
			,panel.grid.major = element_blank()
			,axis.ticks = element_blank()
			,axis.text.y = element_blank()
			,axis.text.x = element_blank()
			,axis.title.y = element_blank()
			,axis.title.x = element_blank()
			,legend.position='none'
	)
	
	size_scale <- scale_size(limits = c(0,1),range=cor_text_limits)
	the.plot<-ggplot(data=z)
	if(type=="bins")
		the.plot<-the.plot+bin_layer
	else if(type=="points")
		the.plot<-the.plot+points_layer + scale_alpha_identity()
	the.plot<-the.plot+var_text+
			cor_text+
			f+
			o+
			size_scale
	if(type=="bins")
		the.plot<-the.plot+scale_fill_gradient(low="grey", high="black")
	if(lines)
		the.plot<-the.plot+lm_line_layer
	if(main=="auto")
		main<-cor.mat[[1]][[1]]$method
	the.plot<-the.plot+ggtitle(main)
	return(the.plot)
}
