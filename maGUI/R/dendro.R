display<-function(h,...){
	tr<-tr
	addHandlerDoubleclick(tr,handler=function(h,...){
		tree_x<-h$obj[]
		if(length(tree_x)>1)
		{
			if(tree_x[1]=="Affymetrix")
			{
				try(({dat2Affy.m<-dat2Affy.m;aqc<-aqc;dat2Affy.f<-dat2Affy.f;dat2Affy.s<-dat2Affy.s;DE_Affy<-DE_Affy;
				pca_Affy<-pca_Affy;sample.dist_Affy<-sample.dist_Affy;Clas_Affy<-Clas_Affy;
				GOresultBP_Affy<-GOresultBP_Affy;GOresultMF_Affy<-GOresultMF_Affy;GOresultCC_Affy<-GOresultCC_Affy;
				KEGGresult_Affy<-KEGGresult_Affy;
				GOtable.outBP_Affy<-GOtable.outBP_Affy;GOtable.outMF_Affy<-GOtable.outMF_Affy;
				GOtable.outCC_Affy<-GOtable.outCC_Affy;KEGGtable.out_Affy<-KEGGtable.out_Affy;
				genes_Affy<-genes_Affy;size_Affy<-size_Affy;
				myGraph_Affy<-myGraph_Affy;}),silent=TRUE)
				
				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(dat2Affy.m)
					disp1<-cbind(Identifier,dat2Affy.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					err<-try(plot(aqc),silent=TRUE)
					if(length(grep("Error",err))!=0){
						plot(dat2Affy.m)
					}
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(dat2Affy.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(dat2Affy.s)
					disp1<-cbind(Identifier,dat2Affy.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(DE_Affy)
					disp1<-cbind(Identifier,DE_Affy)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_Affy,main="Affymetrix PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_Affy,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_Affy,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_Affy))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_Affy,"hypergeoBP_Affymetrix.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_Affy))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_Affy,"hypergeoMF_Affymetrix.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_Affy))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_Affy,"hypergeoCC_Affymetrix.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_Affy))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_Affy,"hypergeoKEGG_Affymetrix.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						Identifier<-rownames(GOtable.outBP_Affy)
						disp1<-cbind(Identifier,GOtable.outBP_Affy)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_Affy)
						disp1<-cbind(Identifier,GOtable.outMF_Affy)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_Affy)
						disp1<-cbind(Identifier,GOtable.outCC_Affy)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(KEGGtable.out_Affy)
					disp1<-cbind(Identifier,KEGGtable.out_Affy)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_Affy,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_Affy,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_Affy,nodeAttrs=makeNodeAttrs(myGraph_Affy,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Agilent_OneColor")
			{
				try(({datAgOne2.m<-datAgOne2.m;datAgOne2.f<-datAgOne2.f;datAgOne2.s<-datAgOne2.s;DE_Ag1<-DE_Ag1;
				pca_Ag1<-pca_Ag1;sample.dist_Ag1<-sample.dist_Ag1;Clas_Ag1<-Clas_Ag1;
				GOresultBP_Ag1<-GOresultBP_Ag1;GOresultMF_Ag1<-GOresultMF_Ag1;GOresultCC_Ag1<-GOresultCC_Ag1;
				KEGGresult_Ag1<-KEGGresult_Ag1;
				GOtable.outBP_Ag1<-GOtable.outBP_Ag1;GOtable.outMF_Ag1<-GOtable.outMF_Ag1;
				GOtable.outCC_Ag1<-GOtable.outCC_Ag1;KEGGtable.out_Ag1<-KEGGtable.out_Ag1;
				genes_Ag1<-genes_Ag1;size_Ag1<-size_Ag1;
				myGraph_Ag1<-myGraph_Ag1;}),silent=TRUE)
				
				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(datAgOne2.m)
					disp1<-cbind(Identifier,datAgOne2.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					boxplot(datAgOne2.m)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(datAgOne2.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(datAgOne2.s)
					disp1<-cbind(Identifier,datAgOne2.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_Ag1)
					disp1<-cbind(Identifier,DE_Ag1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_Ag1,main="Agilent_OneColor PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_Ag1,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_Ag1,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_Ag1))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_Ag1,"hypergeoBP_Agilent_OneColor.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_Ag1))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_Ag1,"hypergeoMF_Agilent_OneColor.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_Ag1))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_Ag1,"hypergeoCC_Agilent_OneColor.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_Ag1))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_Ag1,"hypergeoKEGG_Agilent_OneColor.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_Ag1)
						disp1<-cbind(Identifier,GOtable.outBP_Ag1)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_Ag1)
						disp1<-cbind(Identifier,GOtable.outMF_Ag1)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_Ag1)
						disp1<-cbind(Identifier,GOtable.outCC_Ag1)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_Ag1)
					disp1<-cbind(Identifier,KEGGtable.out_Ag1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_Ag1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_Ag1,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_Ag1,nodeAttrs=makeNodeAttrs(myGraph_Ag1,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Agilent_TwoColor")
			{
				try(({datAgTwo2.m<-datAgTwo2.m;datAgTwo2.f<-datAgTwo2.f;datAgTwo2.s<-datAgTwo2.s;DE_Ag2<-DE_Ag2;
				pca_Ag2<-pca_Ag2;sample.dist_Ag2<-sample.dist_Ag2;Clas_Ag2<-Clas_Ag2;
				GOresultBP_Ag2<-GOresultBP_Ag2;GOresultMF_Ag2<-GOresultMF_Ag2;GOresultCC_Ag2<-GOresultCC_Ag2;
				KEGGresult_Ag2<-KEGGresult_Ag2;
				GOtable.outBP_Ag2<-GOtable.outBP_Ag2;GOtable.outMF_Ag2<-GOtable.outMF_Ag2;
				GOtable.outCC_Ag2<-GOtable.outCC_Ag2;KEGGtable.out_Ag2<-KEGGtable.out_Ag2;
				genes_Ag2<-genes_Ag2;size_Ag2<-size_Ag2;
				myGraph_Ag2<-myGraph_Ag2;}),silent=TRUE)

				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(datAgTwo2.m)
					disp1<-cbind(Identifier,datAgTwo2.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					boxplot(datAgTwo2.m)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(datAgTwo2.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(datAgTwo2.s)
					disp1<-cbind(Identifier,datAgTwo2.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_Ag2)
					disp1<-cbind(Identifier,DE_Ag2)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_Ag2,main="Agilent_TwoColor PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_Ag2,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_Ag2,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_Ag2))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_Ag2,"hypergeoBP_Agilent_TwoColor.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_Ag2))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_Ag2,"hypergeoMF_Agilent_TwoColor.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{					
							if(dim(summary(GOresultCC_Ag2))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_Ag2,"hypergeoCC_Agilent_TwoColor.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_Ag2))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_Ag2,"hypergeoKEGG_Agilent_TwoColor.html",append=T)
						}
					}

				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_Ag2)
						disp1<-cbind(Identifier,GOtable.outBP_Ag2)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_Ag2)
						disp1<-cbind(Identifier,GOtable.outMF_Ag2)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_Ag2)
						disp1<-cbind(Identifier,GOtable.outCC_Ag2)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_Ag2)
					disp1<-cbind(Identifier,KEGGtable.out_Ag2)
					disp<-gtable(KEGGtable.out_Ag2,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_Ag2,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_Ag2,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_Ag2,nodeAttrs=makeNodeAttrs(myGraph_O,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Illumina_Beadarray")
			{
				try(({datIllBA2.m2<-datIllBA2.m2;datIllBA2.f<-datIllBA2.f;datIllBA2.s<-datIllBA2.s;DE_Il_B<-DE_Il_B;
				pca_Il_B<-pca_Il_B;sample.dist_Il_B<-sample.dist_Il_B;Clas_Il_B<-Clas_Il_B;
				GOresultBP_Il_B<-GOresultBP_Il_B;GOresultMF_Il_B<-GOresultMF_Il_B;GOresultCC_Il_B<-GOresultCC_Il_B;
				KEGGresult_Il_B<-KEGGresult_Il_B;
				GOtable.outBP_Il_B<-GOtable.outBP_Il_B;GOtable.outMF_Il_B<-GOtable.outMF_Il_B;
				GOtable.outCC_Il_B<-GOtable.outCC_Il_B;KEGGtable.out_Il_B<-KEGGtable.out_Il_B;
				genes_Il_B<-genes_Il_B;size_Il_B<-size_Il_B;
				myGraph_Il_B<-myGraph_Il_B;}),silent=TRUE)

				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(datIllBA2.m2)
					disp1<-cbind(Identifier,datIllBA2.m2)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					boxplot(datIllBA2.m2)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(datIllBA2.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(datIllBA2.s)
					disp1<-cbind(Identifier,datIllBA2.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_Il_B)
					disp1<-cbind(Identifier,DE_Il_B)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_Il_B,main="Illumina_Beadarray PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_Il_B,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_Il_B,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_Il_B))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_Il_B,"hypergeoBP_Illumina_Beadarray.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_Il_B))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_Il_B,"hypergeoMF_Illumina_Beadarray.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_Il_B))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_Il_B,"hypergeoCC_Illumina_Beadarray.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_Il_B))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_Il_B,"hypergeoKEGG_Illumina_Beadarray.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_Il_B)
						disp1<-cbind(Identifier,GOtable.outBP_Il_B)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_Il_B)
						disp1<-cbind(Identifier,GOtable.outMF_Il_B)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_Il_B)
						disp1<-cbind(Identifier,GOtable.outCC_Il_B)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_Il_B)
					disp1<-cbind(Identifier,KEGGtable.out_Il_B)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_Il_B,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_Il_B,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_Il_B,nodeAttrs=makeNodeAttrs(myGraph_Il_B,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Illumina_Lumi")
			{
				try(({lumi_NQ.m<-lumi_NQ.m;lumi_NQ.f<-lumi_NQ.f;lumi_NQ.s<-lumi_NQ.s;DE_Il_L<-DE_Il_L;
				pca_Il_L<-pca_Il_L;sample.dist_Il_L<-sample.dist_Il_L;Clas_Il_L<-Clas_Il_L;
				GOresultBP_Il_L<-GOresultBP_Il_L;GOresultMF_Il_L<-GOresultMF_Il_L;GOresultCC_Il_L<-GOresultCC_Il_L;
				KEGGresult_Il_L<-KEGGresult_Il_L;
				GOtable.outBP_Il_L<-GOtable.outBP_Il_L;GOtable.outMF_Il_L<-GOtable.outMF_Il_L;
				GOtable.outCC_Il_L<-GOtable.outCC_Il_L;KEGGtable.out_Il_L<-KEGGtable.out_Il_L;
				genes_Il_L<-genes_Il_L;size_Il_L<-size_Il_L;
				myGraph_Il_L<-myGraph_Il_L;}),silent=TRUE)
				
				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(lumi_NQ.m)
					disp1<-cbind(Identifier,lumi_NQ.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(lumi_NQ.m,what="boxplot")
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(lumi_NQ.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(lumi_NQ.s)
					disp1<-cbind(Identifier,lumi_NQ.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_Il_L)
					disp1<-cbind(Identifier,DE_Il_L)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_Il_L,main="Illumina_Lumi PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_Il_L,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_Il_L,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_Il_L))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_Il_L,"hypergeoBP_Illumina_Lumi.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_Il_L))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_Il_L,"hypergeoMF_Illumina_Lumi.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_Il_L))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_Il_L,"hypergeoCC_Illumina_Lumi.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_Il_L))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_Il_L,"hypergeoKEGG_Illumina_Lumi.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_Il_L)
						disp1<-cbind(Identifier,GOtable.outBP_Il_L)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_Il_L)
						disp1<-cbind(Identifier,GOtable.outMF_Il_L)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_Il_L)
						disp1<-cbind(Identifier,GOtable.outCC_Il_L)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_Il_L)
					disp1<-cbind(Identifier,KEGGtable.out_Il_L)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_Il_L,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_Il_L,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_Il_L,nodeAttrs=makeNodeAttrs(myGraph_Il_L,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Nimblegen")
			{
				try(({data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m;
				data.matrix_Nimblegen2.f<-data.matrix_Nimblegen2.f;data.matrix_Nimblegen2.s<-data.matrix_Nimblegen2.s;
				DE_N<-DE_N;pca_N<-pca_N;sample.dist_N<-sample.dist_N;Clas_N<-Clas_N;
				GOresultBP_N<-GOresultBP_N;GOresultMF_N<-GOresultMF_N;GOresultCC_N<-GOresultCC_N;KEGGresult_N<-KEGGresult_N;
				GOtable.outBP_N<-GOtable.outBP_N;GOtable.outMF_N<-GOtable.outMF_N;GOtable.outCC_N<-GOtable.outCC_N;
				KEGGtable.out_N<-KEGGtable.out_N;genes_N<-genes_N;size_N<-size_N;
				myGraph_N<-myGraph_N;}),silent=TRUE)
				
				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(data.matrix_Nimblegen2.m)
					disp1<-cbind(Identifier,data.matrix_Nimblegen2.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					boxplot(data.matrix_Nimblegen2.m)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(data.matrix_Nimblegen2.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(data.matrix_Nimblegen2.s)
					disp1<-cbind(Identifier,data.matrix_Nimblegen2.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_N)
					disp1<-cbind(Identifier,DE_N)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_N,main="Nimblegen PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_N,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_N,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_N))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_N,"hypergeoBP_Nimblegen.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_N))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_N,"hypergeoMF_Nimblegen.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_N))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_N,"hypergeoCC_Nimblegen.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_N))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_N,"hypergeoKEGG_Nimblegen.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_N)
						disp1<-cbind(Identifier,GOtable.outBP_N)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_N)
						disp1<-cbind(Identifier,GOtable.outMF_N)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_N)
						disp1<-cbind(Identifier,GOtable.outCC_N)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_N)
					disp1<-cbind(Identifier,KEGGtable.out_N)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_N,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_N,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_N,nodeAttrs=makeNodeAttrs(myGraph_N,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Series_Matrix")
			{
				try(({data.matrixNorm.m<-data.matrixNorm.m;data.matrixNorm.f<-data.matrixNorm.f;
				data.matrixNorm.s<-data.matrixNorm.s;DE_S<-DE_S;pca_S<-pca_S;sample.dist_S<-sample.dist_S;Clas_S<-Clas_S;
				GOresultBP_S<-GOresultBP_S;GOresultMF_S<-GOresultMF_S;GOresultCC_S<-GOresultCC_S;KEGGresult_S<-KEGGresult_S;
				GOtable.outBP_S<-GOtable.outBP_S;GOtable.outMF_S<-GOtable.outMF_S;GOtable.outCC_S<-GOtable.outCC_S;
				KEGGtable.out_S<-KEGGtable.out_S;genes_S<-genes_S;size_S<-size_S;
				myGraph_S<-myGraph_S;}),silent=TRUE)				
				
				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(data.matrixNorm.m)
					disp1<-cbind(Identifier,data.matrixNorm.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					boxplot(data.matrixNorm.m)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(data.matrixNorm.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(data.matrixNorm.s)
					disp1<-cbind(Identifier,data.matrixNorm.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_S)
					disp1<-cbind(Identifier,DE_S)	
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_S,main="Series_Matrix PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_S,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_S,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_S))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_S,"hypergeoBP_Series_Matrix.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_S))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_S,"hypergeoMF_Series_Matrix.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_S))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_S,"hypergeoCC_Series_Matrix.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_S))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_S,"hypergeoKEGG_Series_Matrix.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_S)
						disp1<-cbind(Identifier,GOtable.outBP_S)	
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_S)
						disp1<-cbind(Identifier,GOtable.outMF_S)	
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_S)
						disp1<-cbind(Identifier,GOtable.outCC_S)	
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_S)
					disp1<-cbind(Identifier,KEGGtable.out_S)	
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_S,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_S,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_S,nodeAttrs=makeNodeAttrs(myGraph_S,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		else if(tree_x[1]=="Online_Data")
			{
				try(({data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m;data.matrix_onlineNorm.f<-data.matrix_onlineNorm.f;
				data.matrix_onlineNorm.s<-data.matrix_onlineNorm.s;DE_O<-DE_O;
				pca_O<-pca_O;sample.dist_O<-sample.dist_O;Clas_O<-Clas_O;
				GOresultBP_O<-GOresultBP_O;GOresultMF_O<-GOresultMF_O;GOresultCC_O<-GOresultCC_O;KEGGresult_O<-KEGGresult_O;
				GOtable.outBP_O<-GOtable.outBP_O;GOtable.outMF_O<-GOtable.outMF_O;GOtable.outCC_O<-GOtable.outCC_O;
				KEGGtable.out_O<-KEGGtable.out_O;genes_O<-genes_O;size_O<-size_O;
				myGraph_O<-myGraph_O;}),silent=TRUE)

				if(tree_x[2]=="Normalization")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(data.matrix_onlineNorm.m)
					disp1<-cbind(Identifier,data.matrix_onlineNorm.m)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="QC_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					boxplot(data.matrix_onlineNorm.m)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Filtered")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					ps1<-as.data.frame(data.matrix_onlineNorm.f)	
					Identifier<-rownames(ps1)
					disp1<-cbind(Identifier,ps1)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Stat_Significant")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					Identifier<-rownames(data.matrix_onlineNorm.s)
					disp1<-cbind(Identifier,data.matrix_onlineNorm.s)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="DGE")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(DE_O)
					disp1<-cbind(Identifier,DE_O)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="PCA_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					try(plot(pca_O,main="Online_Data PCA"),silent=TRUE)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Cluster_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					plot(hclust(sample.dist_O,method="complete"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Classification")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE	
					heatmap(Clas_O,Rowv=NA)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="GSEA_GO")
				{
					if(tree_x[3]=="BP")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO BP",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultBP_O))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO BP",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO BP",delay=5,parent=c(600,400))
								htmlReport(GOresultBP_O,"hypergeoBP_Online_Data.html",append=T)
							}
						}
					}
					if(tree_x[3]=="MF")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO MF",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultMF_O))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO MF",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO MF",delay=5,parent=c(600,400))
								htmlReport(GOresultMF_O,"hypergeoMF_Online_Data.html",append=T)
							}
						}
					}
					if(tree_x[3]=="CC")
					{
						html1<-gconfirm("Generate Html Report?","GSEA GO CC",icon="question")
						if(html1==TRUE)
						{
							if(dim(summary(GOresultCC_O))[1]==0)
							{
								galert("No reports to generate",title="GSEA GO CC",delay=5,parent=c(600,400))
							}
							else
							{
								galert("Generating Html Report","GSEA GO CC",delay=5,parent=c(600,400))
								htmlReport(GOresultCC_O,"hypergeoCC_Online_Data.html",append=T)
							}
						}
					}
				}
				if(tree_x[2]=="GSEA_KEGG")
				{
					html1<-gconfirm("Generate Html Report?","GSEA KEGG",icon="question")
					if(html1==TRUE)
					{
						if(dim(summary(KEGGresult_O))[1]==0)
						{
							galert("No reports to generate",title="GSEA KEGG",delay=5,parent=c(600,400))
						}
						else
						{
							galert("Generating Html Report",title="GSEA KEGG",delay=5,parent=c(600,400))
							htmlReport(KEGGresult_O,"hypergeoKEGG_Online_Data.html",append=T)
						}
					}
				}
				if(tree_x[2]=="GSTA_GO")
				{
					if(tree_x[3]=="BP")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outBP_O)
						disp1<-cbind(Identifier,GOtable.outBP_O)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="MF")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outMF_O)
						disp1<-cbind(Identifier,GOtable.outMF_O)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
					if(tree_x[3]=="CC")
					{
						galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
						delete(g2,g2_1)
						g2_1<<-ggroup(container=g2,horizontal=FALSE)
						visible(g2)<-FALSE	
						Identifier<-rownames(GOtable.outCC_O)
						disp1<-cbind(Identifier,GOtable.outCC_O)
						disp<-gtable(disp1,container=g2_1)
						size(disp)<-c(650,450)
						visible(g2)<-TRUE
						galert("Done",title="Loading",delay=5,parent=c(600,400))
					}
				}
				if(tree_x[2]=="GSTA_KEGG")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE	
					Identifier<-rownames(KEGGtable.out_O)
					disp1<-cbind(Identifier,KEGGtable.out_O)
					disp<-gtable(disp1,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Identifier_Symbol")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					visible(g2)<-FALSE
					disp<-gtable(genes_O,container=g2_1)
					size(disp)<-c(650,450)
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="SSE_Plot")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					ssize.plot(size_O,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
				if(tree_x[2]=="Coexpression_Network")
				{
					galert("Plz.. wait while loading",title="Loading",delay=5,parent=c(600,400))
					delete(g2,g2_1)
					g2_1<<-ggroup(container=g2,horizontal=FALSE)
					size(g2_1)<-c(650,400)
					plotarea=NULL;
					plotarea<<-ggraphics(ps=5,use.scrollwindow=TRUE,horizontal=FALSE,container=g2_1)
					visible(g2)<-FALSE
					disp<-plot(myGraph_O,nodeAttrs=makeNodeAttrs(myGraph_O,fontsize=18,fillcolor="grey"))
					visible(g2)<-TRUE
					galert("Done",title="Loading",delay=5,parent=c(600,400))
				}
			}
		}
	}	
	)
}
