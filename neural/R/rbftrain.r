rbftrain<-function(inp,neurons,out,weight=c(),dist=c(),alfa=0.2,it=40,err=0,sigma=NaN,online=TRUE,permute=TRUE,visual=TRUE, ...){

		rect<-function(x,y,value){
			xx<-c(x,x+30,x+30,x)
			yy<-c(y,y,y+30,y+30)
			polygon(xx,yy,col="green");
		}

		drawnet<-function(run=FALSE){
			plot(1:700,1:700,xlab="",ylab="",type="n",axes=FALSE);
			title(main="RBF network");

			polygon(c(0,65,65,0),c(125,125,87,87),col="lightblue2");
			text(32,109,"Prev",cex=0.9);
			polygon(c(80,145,145,80),c(125,125,87,87),col="lightblue2");
			text(112,109,"Next",cex=0.9);
			polygon(c(640,720,720,640),c(650,650,700,700),col="lightblue2");
			text(680,678,ifelse(run==FALSE,"START","EXIT"),cex=0.9);
			if ((run)&(ncol(inp)==1)&(ncol(out)==1)){
				polygon(c(640,720,720,640),c(550,550,600,600),col="lightblue2");
				text(680,578,"DRAW",cex=0.9);
			}
			text(25,70,"Iteration:",cex=0.9);
			polygon(c(75,160,160,75),c(45,45,85,85),col="turquoise");
			polygon(c(160,240,240,160),c(45,45,85,85),col="snow3");
			text(170,67,ifelse(run==FALSE,0,iter),pos=2,cex=0.9);
			text(245,67,it,pos=2,cex=0.9);

			text(25,30,"Error:",cex=0.9);
			polygon(c(75,160,160,75),c(5,5,45,45),col="turquoise");
			polygon(c(160,240,240,160),c(5,5,45,45),col="snow3");
			text(170,27,round(error*10^(4-max(round(log10(error)),0)))/10^(4-max(round(log10(error)),0)),pos=2,cex=0.9);
			text(245,27,err,pos=2,cex=0.9);


			for(i in 1:length(cordx)){
				for(j in 1:length(cordx[[i]])){
					rect(cordx[[i]][j],cordy[[i]][j],val[[i]][j]);
					if (i==1)
						text(cordx[[i]][j]-10,cordy[[i]][j]+5,round(val[[i]][j]*100)/100,cex=0.8)
					else text(cordx[[i]][j]+40,cordy[[i]][j]+5,round(val[[i]][j]*100)/100,cex=0.8);
				}
			}
			for(i in 1:(length(cordx)-1)){
				for(j in 1:length(cordx[[i]])){
					for(k in 1:length(cordx[[i+1]])){
						xx<-c(cordx[[i]][j]+30,cordx[[i+1]][k]);
						yy<-c(cordy[[i]][j]+15,cordy[[i+1]][k]+15);
						lines(xx,yy,col="blue");
						text((xx[1]+xx[2])/2,(yy[1]+yy[2])/2,round(weight[[i]][j,k]*100)/100,cex=0.8)
					}
				}
			}
		}


		graph<-function(){
			plot(1:700,1:700,xlab="",ylab="",type="n",axes=FALSE);

			dx<-(max(inp)-min(inp))*1.2
			dy<-(max(out)-min(out))*1.2
			minx<-min(inp)-dx/12
			miny<-min(out)-dy/12

			title(main="RBF network");
			polygon(c(640,720,720,640),c(650,650,700,700),col="lightblue2");
			text(680,678,"EXIT",cex=0.9);
			polygon(c(640,720,720,640),c(550,550,600,600),col="lightblue2");
			text(680,578,"NETW",cex=0.9);

			lines(c(40,40),c(20,720));
			lines(c(40,740),c(20,20));

			for(i in 0:7){
				lines(c(30,50),c(100*i+20,100*i+20));
				text(x=0,y=i*100+20,round(i/7*dy+miny,digits=3),cex=0.7);
				lines(c(100*i+40,100*i+40),c(10,30));
				text(x=i*100+40,y=0,round(i/7*dx+minx,digits=3),cex=0.7);
			}


			for(i in 0:700) lines(c(i+40,i+41),c((valuate3(minx+i/700*dx)-miny)/dy*700+20,(valuate3(minx+(i+1)/700*dx)-miny)/dy*700+20))
			for(i in 1:length(inp)) points((inp[i]-minx)/dx*700+40,(valuate3(inp[i])-miny)/dy*700+20,col="red2")
			for(i in 1:length(inp)) points((inp[i]-minx)/dx*700+40,(out[i]-miny)/dy*700+20,col="purple1")

		}

		gauss<-function(x,sigma) exp(-(x^2)/(2*sigma^2));

		ident<-function(x) x;

		identdif<-function(x) 1;
	
		v<-function(l,j) {vv<-0;for(i in 1:ls[l-1]) vv<-vv+val[[l-1]][i]*weight[[l-1]][i,j];vv+dist[j]};


		valuate1<-function(watch){
			value<-vector(1,mode="list");
			value[[1]]<-inp[watch,];
			total<-0;
			ee<-c();


			for(j in 1:ls[2]){
				e<-0;
				for(k in 1:ls[1]) 
					e<-e+abs(value[[1]][k]-weight[[1]][k,j]);
				total<-total+gauss(e,sigmavalue[k,j]);
			}

			if ((total==0)||(is.nan(total))){
				ee<-rep(0,ls[2]);
			}
			else{
				for(j in 1:ls[2]){
					e<-0;
					for(k in 1:ls[1])
						e<-e+abs(value[[1]][k]-weight[[1]][k,j]);
					ee<-c(ee,gauss(e,sigmavalue[k,j])/total);
				}
			}
			ee;
		}

		valuate2<-function(centre){
			value<-vector(3,mode="list");
			value[[1]]<-inp[watch,];
			value[[2]]<-centre[watch,];
			
			ee<-c();
			for(j in 1:ls[3]){
				e<-0;
				for(k in 1:ls[2]) 
					e<-e+value[[2]][k]*weight[[2]][k,j];
				ee<-c(ee,ident(e+dist[j]));
			}
			value[[3]]<-ee;
			
			value;
		}

		valuate3<-function(pont){
			value<-vector(2,mode="list");
			value[[1]]<-pont;

			total<-0;
			ee<-c();


			for(j in 1:ls[2]){
				e<-0;
				for(k in 1:ls[1]) 
					e<-e+abs(value[[1]][k]-weight[[1]][k,j]);
				total<-total+gauss(e,sigmavalue[k,j]);
			}


			if ((total==0)||(is.nan(total))){
				ee<-rep(0,ls[2])
			}
			else{
				for(j in 1:ls[2]){
					e<-0;
					for(k in 1:ls[1]) 
						e<-e+abs(value[[1]][k]-weight[[1]][k,j]);
					ee<-c(ee,gauss(e,sigmavalue[k,j])/total);
				}
			}
			value[[2]]<-ee;

			e<-0
			for(k in 1:ls[2]) 
				e<-e+value[[2]][k]*weight[[2]][k,1];
			ident(e+dist[1]);
		}


		deltaz<-function(){
			deltak<-vector(3,mode="list");
				dd<-c()
				for(j in 1:ls[3]){
					dd<-c(dd,(out[watch,j]-val[[3]][j]));#*identdif(v(i,j)));
				}
				deltak[[3]]<-dd;
			deltak;
		}


		clust<-function(sampl,k,iter=20){
	
			kozep<-rep(sampl[1],k)
			set<-rep(1,length(sampl))
		
			for(i in 1:k)
				kozep[i]<-(max(sampl)-min(sampl))/k*i+min(sampl)

			valt<-1;szamol<-0;
			while((valt!=0)&(szamol<iter)){
				valt<-0;szamol<-szamol+1;
				for(i in 1:length(sampl)){
					max<-abs(sampl[i]-kozep[set[i]]);
					for(j in 1:k){
						if (abs(sampl[i]-kozep[j])<max){
							valt<-1;
							max<-abs(sampl[i]-kozep[j]);
							set[i]<-j;
						}
					}
				}
				for(i in 1:k){
					summ<-0;cnt<-0;
					for(j in 1:length(sampl))
						if (set[j]==i){
							summ<-summ+sampl[j];cnt<-cnt+1;
						}
					kozep[i]<-ifelse(cnt==0,0,summ/cnt);
				}
			}

			kozep;
		}



		sigmaz<-function(){
				sigmas<-matrix(0,ls[1],ls[2]);
				for(k in 1:ls[1]){

					for(i in 1:ls[2]){
						d<-0
						for(j in 1:ls[2])
							if (i!=j)
								d<-c(d,weight[[1]][k,j]-weight[[1]][k,i]);
						maxim<-max(d);
						minim<-min(d);
						for(j in 1:length(d)){
							if ((d[j]<maxim)&(d[j]>0)) maxim<-d[j]
							if ((d[j]>minim)&(d[j]<0)) minim<-d[j]
						}
						s<-max(abs(minim),maxim)/2		
						sigmas[k,i]<-s*1.1
					}
				}
				sigmas;
		}


	
	if ((length(neurons)!=1)|(neurons<=0)|(neurons%%1!=0)) return("Neuron must be one positive integer number");
	err<-abs(err)
	ls<-c(ncol(inp),neurons,ncol(out))
	if (nrow(inp)!=nrow(out)) return("Different input and output sample length");
	retr<-FALSE;
	retr2<-FALSE;

	if (length(weight)!=0){
			if (length(weight)!=2) return("The weight arguments length must be 2.")
			for(i in 1:2){
				if (nrow(weight[[i]])!=ls[i]) return(paste("The number of rows is different in weight[",i,"].",sep=""));
				if (ncol(weight[[i]])!=ls[i+1]) return(paste("The number of column is different in weight[",i,"].",sep=""))
			}
			retr<-TRUE
		}
	else{
		weight<-vector(2,mode="list")
		weight[[1]]<-matrix(c(0),ls[1],ls[2])
		weight[[2]]<-matrix(runif(ls[2]*ls[3],min=-1),ls[2],ls[3])
	}

	if (length(dist)!=0){
		if (length(dist)!=ls[3]) return("The length of the distortion and the number of neurons in the third layer is different")
		retr2<-TRUE
	}
	else{
		dist<-rep(1,ls[3])
		for(j in 1:ls[3]) dist[j]<-runif(1,min=-1,max=1);
	}
        sigmavalue<-matrix(0,ls[1],ls[2])

	if (visual){
		cordx<-vector(3,mode="list");cordy<-vector(3,mode="list");
		for(i in 1:length(ls)){
			xc<-c();yc<-c();
			for(j in 1:ls[i]){
				xc<-c(xc,350-length(ls)*80+i*160);
				yc<-c(yc,350+ls[i]*30-j*60);
			}
			cordx[[i]]<-xc;cordy[[i]]<-yc;
		}
	}

	

	centre<-matrix(1,nrow(inp),ls[2])
	for(watch in 1:nrow(inp))	centre[watch,]<-valuate1(watch);

	watch<-1;
	val<-valuate2(centre);
	run<-FALSE;
	error<-NaN;
	if (visual) drawnet(run);
	
	ext<-FALSE;valt<-FALSE;graf<-FALSE;
	while(!ext){		
		if (visual){ coor<-locator(1);
			if ((coor$x>640)&(coor$x<720)&(coor$y>550)&(coor$y<620)&(run)&(ncol(inp)==1)&(ncol(out)==1)){
				graf<-!graf;
				if (graf) graph()
				   else drawnet(run);
			}
		}
		else{
			coor<-list(x=0,y=0);
		}

		if ((!visual)|(coor$x>640)&(coor$x<720)&(coor$y>650)&(coor$y<720)){
			if (!run) {
				if (!retr){
					if (visual){
						polygon(c(150,550,550,150),c(400,400,300,300),col="ivory2");
						text(350,350,"Clustering",cex=2);
					}

					if (nrow(inp)<=neurons){
						for(i in 1:ncol(inp))
							weight[[1]][i,]<-clust(inp[,i],neurons);
					}
					else
						for(i in 1:ncol(inp))
							weight[[1]][i,]<-kmeans(inp[,i],neurons)$centers[,1];

					}				
				
				if (is.matrix(sigma)) sigmavalue<-sigma
					else if (is.nan(sigma)) sigmavalue<-sigmaz()
						else if (length(sigma)==1) sigmavalue<-matrix(sigma,ls[1],ls[2])
							else return("The sigma argument is incorrect!!!");
						    
				
				for(i in 1:nrow(inp)) centre[i,]<-valuate1(i);
				val<-valuate2(centre);
				if (visual) drawnet(run);
				iter=0;				

				if (!is.na(it)&(it!=0))
				while ((iter<it)&((is.na(error))|(error>err))){
					iter<-iter+1;
					if (visual){
						polygon(c(75,160,160,75),c(45,45,85,85),col="turquoise");
						text(170,67,iter,pos=2,cex=0.9);
						polygon(c(75,160,160,75),c(5,5,45,45),col="turquoise");
						text(170,27,round(error*10^(4-max(round(log10(error)),0)))/10^(4-max(round(log10(error)),0)),pos=2,cex=0.9);
					}
			
					if (permute) perm<-sample(nrow(inp),nrow(inp)) else perm<-1:nrow(inp);
					w2<-weight;
					t2<-dist;
					error<-0
					for(ii in 1:nrow(inp)){
						watch<-perm[ii];
						val<-valuate2(centre);
						for(j in 1:ncol(out)) error<-error+abs(ifelse(is.na(val[[3]][j]),0,val[[3]][j])-out[watch,j]);
						delta<-deltaz();
						for(k in (length(ls)-1):2){
							for(i in 1:ls[k])
								for(j in 1:ls[k+1]){
									valtoz<-alfa*delta[[k+1]][j]*val[[k]][i];
									w2[[k]][i,j]<-w2[[k]][i,j]+valtoz;
								}							
							for(j in 1:ls[k+1])
								if (online) dist[j]<-dist[j]+alfa*delta[[k+1]][j]
								else t2[[k+1]][j]<-t2[[k+1]][j]+alfa*delta[[k+1]][j];
						}
						if (online) weight<-w2;
					}
					error<-error/(nrow(inp)*ncol(out));
					if (!online){ weight<-w2;dist<-t2;}
				}

		
				watch<-1;
				run<-TRUE;val<-valuate2(centre);
				if (visual) drawnet(run)
				else ext<-TRUE;
			}
			else ext<-TRUE;
		}

		if (visual){
			if ((coor$x>0)&(coor$x<65)&(coor$y>87)&(coor$y<125)&(!graf)){
				watch<-ifelse(watch>1,watch-1,1);val<-valuate2(centre);drawnet(run);
			}
			if ((coor$x>80)&(coor$x<145)&(coor$y>87)&(coor$y<125)&(!graf)){
				watch<-ifelse(watch<nrow(inp),watch+1,nrow(inp));val<-valuate2(centre);drawnet(run);
			}
			if ((valt)&(!graf)) {
				cordx[[vi]][vj]<-coor$x;cordy[[vi]][vj]<-coor$y;valt=FALSE;
				drawnet(run);
			}
			else
			{
			for(i in 1:length(cordx))
				for(j in 1:length(cordx[[i]]))
					if ((coor$x>cordx[[i]][j])&(coor$x<cordx[[i]][j]+30)&(coor$y>cordy[[i]][j])&(coor$y<cordy[[i]][j]+30)&(!graf)){
						valt<-TRUE;
						vi<-i;vj<-j;
					}
			}
		}
	}

	reslt<-list(weight=weight,dist=dist,neurons=ls,sigma=sigmavalue);
	reslt;
}
