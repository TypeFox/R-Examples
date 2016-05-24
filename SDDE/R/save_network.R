# function save_network
# Save a picture of the augmented network Y with original node in white and  
# additional node in orange.

# Input:
# g1:          network X
# g2:          augmented network Y
# filename:    output filename 

# Output: 
# Two png files. The first one (identified as 1) is the representation of the original network X.
#                The second one (identified as 2) is the representation of the augmented network Y. 

# version 1.0
# created Etienne Lord
# since December 2014
# Layout: layout.fruchterman.reingold, layout.kamada.kawai
# Color from Kelly KL. Twenty-two colors of maximum contrast. Col Eng 1976;3:26-27. See also :
# http://jfly.iam.u-tokyo.ac.jp/color/, and Campadelli, P., Posenato, R., & Schettini, R. (1999). An algorithm for the selection of high-contrast color sets. Color Research & Application, 24(2), 132-138.
	  
save_network<-function(g1,g2, filename, layout=layout.kamada.kawai, taxnames='', mode='png', imagesize=800) {	
	#g1 = blue
	#g2 = orange
	if (mode!="screen") cat("Saving graph [",filename,"]\n");	
	# Color from Kelly KL. Twenty-two colors of maximum contrast. Col Eng 1976;3:26-27. See also :
	# http://jfly.iam.u-tokyo.ac.jp/color/, and Campadelli, P., Posenato, R., & Schettini, R. (1999). An algorithm for the selection of high-contrast color sets. Color Research & Application, 24(2), 132-138.
	Kelly_color=c("#FFFFFF","#FFB300","#803E75","#FF6800","#A6BDD7","#C10020","#CEA262","#817066","#007D34","#F6768E","#00538A","#FF7A5C","#53377A","#FF8E00","#B32851","#F4C800","#7F180D","#93AA00","#593315","#F13A13","#232C16","#000000")	  
	nv=c("Original nodes","Augmented nodes"); #name vector
	cv=c("white","orange"); #color vector
	if (taxnames=='') {
		V(g2)$color <- ifelse((V(g2)$name %in% V(g1)$name), "white", "orange");
		V(g1)$color <-"white";
	} else if (taxnames=='allgroup'&&!is.null(V(g2)$name)&&!is.null(V(g2)$tax)) {					
			i=1;
			nv=c();
			cv=c();
			
			cl=Kelly_color;
			if (length(table(V(g2)$tax))<10) cl=rainbow(length(table(V(g2)$tax))); 
			#Create a table of the coloring
			for (a in names(table(V(g2)$tax))) {			
				cv=c(cv,cl[i]);
				nv=c(nv,a);
				i=i+1;
			}			
			for (i in 1:length(V(g2)$name)) {				
				color_pos=match(V(g2)[i]$tax, nv);								
				V(g2)[i]$color=cv[color_pos];				
			}
			for (i in 1:length(V(g1)$name)) {				
				color_pos=match(V(g1)[i]$tax, nv);								
				V(g1)[i]$color=cv[color_pos];				
			}
	} else if (taxnames!=''&&taxnames!='allgroup'){			
		V(g2)$color <- ifelse((V(g2)$tax != taxnames||V(g2)$name!=taxnames), "white", "orange")
		V(g1)$color <-"white";	
	} else {
		
		V(g2)$color <- ifelse((V(g2)$name %in% V(g1)$name), "white", "orange");
		V(g1)$color <-"white";
	}
	V(g2)$size=10;
	V(g1)$size=10;
	
	#Unique
	vertex_of_g2<-V(g2)[!(V(g2) %in% V(g1))]
	vertex_of_g1=length(V(g1)$name)
	if (length(V(g2)) > 50) {
		save_network_big(g1,g2, filename, layout, taxnames, mode, imagesize);
	} else {		
		if (mode=='png') {
			f2<-paste(filename,"_2",".png", sep="");
			png(f2, imagesize, imagesize, pointsize = 14);
			plot(g2, 
			layout=layout, 
			main=filename,	#specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black', 		#the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=V(g2)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();		
		} else if (mode=='svg') {
			f2<-paste(filename,"_2",".svg", sep="");
			svg(f2);
			plot(g2, 
			layout=layout, 
			main=filename,	#specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black', 		#the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=V(g2)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();		
		} else if (mode=='eps'){
			f2<-paste(filename,"_2",".eps", sep="");
			postscript(f2, horizontal = FALSE, fonts=c("serif", "Palatino"), onefile = FALSE, paper = "special", height = 8, width = 10)
				plot(g2, 
				layout=layout, 
				vertex.label.family="serif",
				edge.label.family="Palatino", 			
				main=filename,	#specifies the title
				vertex.label.dist=0,			#puts the name labels slightly off the dots
				vertex.frame.color='black', 		#the color of the border of the dots 
				vertex.label.color='black',		#the color of the name labels
				vertex.label.font=1,			#the font of the name labels
				vertex.label=V(g1)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
				a=dev.off();		
		} else {			
				plot(g2, 
				layout=layout, 
				vertex.label.family="serif",
				edge.label.family="Palatino", 			
				main=filename,	#specifies the title
				vertex.label.dist=0,			#puts the name labels slightly off the dots
				vertex.frame.color='black', 		#the color of the border of the dots 
				vertex.label.color='black',		#the color of the name labels
				vertex.label.font=1,			#the font of the name labels
				vertex.label=V(g1)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);				
		}		
		#layout=layout,	# the layout method. see the igraph documentation for details
		
		if (mode=='png') {
			f1<-paste(filename,"_1",".png", sep="");
			png(f1, 800, 800, pointsize = 14);
			plot(g1, 
			layout=layout, 
			main=filename,	#specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black', 		#the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=V(g1)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();
		} else if (mode=='svg') {
			f1<-paste(filename,"_1",".svg", sep="");
			svg(f1);
			plot(g1, 
			layout=layout, 
			main=filename,	#specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black', 		#the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=V(g1)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();
		} else if (mode=='eps') {
				f1<-paste(filename,"_1",".eps", sep="");
				postscript(f1, fonts=c("serif", "Palatino"),horizontal = FALSE, onefile = FALSE, paper = "special", height = 8, width = 10)
				plot(g1, 
					layout=layout, 
					vertex.label.family="serif",
					edge.label.family="Palatino", 			
					main=filename,	#specifies the title
					vertex.label.dist=0,			#puts the name labels slightly off the dots
					vertex.frame.color='black', 		#the color of the border of the dots 
					vertex.label.color='black',		#the color of the name labels
					vertex.label.font=1,			#the font of the name labels
					vertex.label=V(g1)$name		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
				a=dev.off();
		} 
		if (mode!="screen") print("done")
	}
}

save_network_big<-function(g1,g2, filename, layout=layout.kamada.kawai, taxnames='', mode='png', imagesize=2500) {		
	#g1 = white
	#g2 = orange
	
	
		Kelly_color=c("#FFFFFF","#FFB300","#803E75","#FF6800","#A6BDD7","#C10020","#CEA262","#817066","#007D34","#F6768E","#00538A","#FF7A5C","#53377A","#FF8E00","#B32851","#F4C800","#7F180D","#93AA00","#593315","#F13A13","#232C16","#000000")	  
	nv=c("Original nodes","Augmented nodes"); #name vector
	cv=c("white","orange"); #color vector
	if (taxnames=='') {
		V(g2)$color <- ifelse((V(g2)$name %in% V(g1)$name), "white", "orange")
	} else if (taxnames=='allgroup'&&!is.null(V(g2)$name)) {		
			i=1;
			nv=c();
			cv=c();
			cl=Kelly_color;
			if (length(table(V(g2)$tax))<10) {
				cl=rainbow(length(table(V(g2)$tax))); 
				# if (length(table(V(g2)$tax))==3) {
					# cl[1]="white";
					# cl[2]="orange";
					# cl[3]="black";
				# }
			}
			#Create a table of the coloring
			for (a in names(table(V(g2)$tax))) {			
				cv=c(cv,cl[i]);
				nv=c(nv,a);
				i=i+1;
			}			
			for (i in 1:length(V(g2)$name)) {				
				color_pos=match(V(g2)[i]$tax, nv);								
				V(g2)[i]$color=cv[color_pos];				
			}
			for (i in 1:length(V(g1)$name)) {				
				color_pos=match(V(g1)[i]$tax, nv);								
				V(g1)[i]$color=cv[color_pos];				
			}
	} else if (taxnames!=''&&taxnames!='allgroup'){		
		V(g2)$color <- ifelse((V(g2)$tax != taxnames||V(g2)$name!=taxnames), "white", "orange")
		V(g1)$color <-"white";	
	} else {
		V(g2)$color <- ifelse((V(g2)$name %in% V(g1)$name), "white", "orange");
		V(g1)$color <-"white";
	}
	V(g2)$size=2.5;
	V(g1)$size=2.5;
	
	#filename2<-paste("graph_randomy",count,".txt", sep="");
	#sink(filename2)   
	#print(res1);
	#print(res2);
	#sink()
	#layout=layout;	# the layout method. see the igraph documentation for details
	if (mode=='png') {
			f2<-paste(filename,"_2",".png", sep="");
			png(f2, imagesize, imagesize, pointsize = 14);
			plot(g2, 
			layout=layout, 
			main=filename,	#specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black', 		#the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();
	} else if (mode=='svg') {
			f2<-paste(filename,"_2",".svg", sep="");
			svg(f2);		
			plot(g2, 
			layout=layout, 
			main=filename,	#specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black', 		#the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();
	} else if (mode=='eps') {
			f2<-paste(filename,"_2",".eps", sep="");
			postscript(f2, , fonts=c("serif", "Palatino"), horizontal = FALSE, onefile = FALSE, paper = "special", height = 8, width = 10)
			plot(g2, 
				layout=layout, 
				vertex.label.family="serif",
				edge.label.family="Palatino", 	
				vertex.label.cex=1,
				main=filename,	#specifies the title
				vertex.label.dist=0,			#puts the name labels slightly off the dots
				vertex.frame.color='black', 		#the color of the border of the dots 
				vertex.label.color='black',		#the color of the name labels
				vertex.label.font=1,			#the font of the name labels
				vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
				a=dev.off();
	} else {
		plot(g2, 
				layout=layout, 
				vertex.label.family="serif",
				edge.label.family="Palatino", 	
				vertex.label.cex=1,
				main=filename,	#specifies the title
				vertex.label.dist=0,			#puts the name labels slightly off the dots
				vertex.frame.color='black', 		#the color of the border of the dots 
				vertex.label.color='black',		#the color of the name labels
				vertex.label.font=1,			#the font of the name labels
				vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
	}
	if (mode=='png') {
			f1<-paste(filename,"_1",".png", sep="");
			png(f1, imagesize, imagesize, pointsize = 14);
			plot(g1, 
			layout=layout, 
			main=filename,	                #specifies the title
			vertex.label.dist=0,			#puts the name labels slightly off the dots
			vertex.frame.color='black',     #the color of the border of the dots 
			vertex.label.color='black',		#the color of the name labels
			vertex.label.font=1,			#the font of the name labels
			vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used			
			);
			legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
			a=dev.off();
	} else if (mode=='svg') {
			f1<-paste(filename,"_1",".svg", sep="");
			svg(f1);		
			plot(g1, 
				layout=layout, 
				vertex.label.family="serif",
				edge.label.family="Palatino", 	
				vertex.label.cex=1,
				main=filename,	#specifies the title
				vertex.label.dist=0,			#puts the name labels slightly off the dots
				vertex.frame.color='black', 		#the color of the border of the dots 
				vertex.label.color='black',		#the color of the name labels
				vertex.label.font=1,			#the font of the name labels
				vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
				a=dev.off();
	} else if (mode=='eps'){
			f1<-paste(filename,"_1",".eps", sep="");
			postscript(f1, , fonts=c("serif", "Palatino"), horizontal = FALSE, onefile = FALSE, paper = "special", height = 8, width = 10)
			plot(g1, 
				layout=layout, 
				vertex.label.family="serif",
				edge.label.family="Palatino", 	
				vertex.label.cex=1,
				main=filename,	#specifies the title
				vertex.label.dist=0,			#puts the name labels slightly off the dots
				vertex.frame.color='black', 		#the color of the border of the dots 
				vertex.label.color='black',		#the color of the name labels
				vertex.label.font=1,			#the font of the name labels
				vertex.label=NA		#specifies the lables of the vertices. in this case the 'name' attribute is used
				);
				legend(x="bottomright",title="Legend", legend=nv,pt.bg=cv, pt.cex=2, col="black",pch=21, yjust=0, lty=0);
				a=dev.off();
	} 
	if (mode!="screen") print("done")
}

#Wrapper to plot to std output
plot_network<-function(g1,g2,layout=layout.kamada.kawai, taxnames='') {	
	save_network(g1,g2,filename='dummy',taxnames=taxnames, layout=layout,mode='screen');
}