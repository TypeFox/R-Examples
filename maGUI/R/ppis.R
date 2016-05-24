ppis<-function(h,...){
ppis_result=NULL;ppis_1=NULL;ppis_2=NULL;f1=NULL;f2=NULL;f_norm1=NULL;f_norm2=NULL;
	w<-gwindow("PPI_Prediction",horizontal=FALSE,height=120,width=250)
	g1<-ggroup(container=w,horizontal=FALSE)
	b1<-gbutton("Normalized_data1",container=g1)
	b2<-gbutton("Normalized_data2",container=g1)
	g2<-ggroup(container=w)
	addHandlerClicked(b1,handler=function(h,...){
		choose_file()
		f1<<-filechoose
		filechoose=NULL
	}
	)
	addHandlerClicked(b2,handler=function(h,...){
		choose_file()
		f2<<-filechoose
		filechoose=NULL
	}
	)
	cancel<-gbutton("CANCEL",border=TRUE,handler=function(h,...){
	dispose(w)
	},container=g2,anchor=c(1,-1))
	ok<-gbutton("OK",border=TRUE,handler=function(h,...){
		dispose(w);
		if(f1=="")galert("Normalized data 1 not selected")
		if(f2=="")galert("Normalized data 2 not selected")
		if(f1=="" && f2=="")galert("Normalized data 1 and 2 not selected")
		if(f1!="" && f2!="")
		{
			f_norm1<<-read.table(f1,sep="\t",row.names=1)
			di1<-dim(f_norm1)
			if(di1[1]>5000){
				if(gconfirm("Higher Dimensions..!\nConsidering only first 5000 rows")==TRUE)
				f_norm1<<-f_norm1[1:50,]
			}
			f_norm2<<-read.table(f2,sep="\t",row.names=1)
			di2<-dim(f_norm2)
			if(di2[1]>5000){
				if(gconfirm("Higher Dimensions..!\nConsidering only first 5000 rows")==TRUE)
				f_norm2<<-f_norm2[1:50,]
			}
			galert("Please wait while PPIs prediction",title="Miscellaneous",delay=5,parent=c(600,400))
			x1<-t(f_norm1)
			z1<-cor(x1,x1,use="everything",method="pearson")
			a<-which(z1[,]>=0.8)
			b<-sqrt(length(z1))
			e=0
			for(i in 1:b)
			{
				l<-e+i
				m<-(e+1):(l-1)
				if(i!=1)
				{
					o<-a[which(a%in%m)]
					p<-o%%b
					if(length(rownames(z1)[p])!=0)
					{
						pp1<-paste(rownames(z1)[i],rownames(z1)[p],sep="-")
						ppis_1<<-c(ppis_1,pp1)
					}
				}
				e<-e+b
			}
			x2<-t(f_norm2)
			z2<-cor(x2,x2,use="everything",method="pearson")
			a<-which(z2[,]>=0.8)
			b<-sqrt(length(z2))
			e=0
			for(i in 1:b)
			{
				l<-e+i
				m<-(e+1):(l-1)
				if(i!=1)
				{
					o<-a[which(a%in%m)]
					p<-o%%b
					if(length(rownames(z2)[p])!=0)
					{
						pp2<-paste(rownames(z2)[i],rownames(z2)[p],sep="-")
						ppis_2<<-c(ppis_2,pp2)
					}
				}
				e<-e+b
			}
			ppis_result<<-ppis_1[which(ppis_1%in%ppis_2)]
			if(length(ppis_result)!=0)
			{
				en<-ginput("Give a name for the results")
				el<-gfile("Select the location to save",type="selectdir")
				setwd(el)
				write.table(ppis_result,file=en,sep="\t")
			}
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
		}
	},container=g2,anchor=c(1,-1))
}
