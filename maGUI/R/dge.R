dge<-function(h,...){
	f<-function(h,...){
		x<<-svalue(h$obj)
		}
		try(({
		dat2Affy.f<-dat2Affy.f;datAgOne2.f<-datAgOne2.f;datAgTwo2.f<-datAgTwo2.f;
		datIllBA2.f<-datIllBA2.f;lumi_NQ.f<-lumi_NQ.f;data.matrix_Nimblegen2.f<-data.matrix_Nimblegen2.f;
		data.matrixNorm.f<-data.matrixNorm.f;data.matrix_onlineNorm.f<-data.matrix_onlineNorm.f;l<-l;tree<-tree;
			}),silent=TRUE)
	platforms=NULL;
	aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;gg=0;hh=0;
	try(({
		if(exists("dat2Affy.f"))aa=length(dat2Affy.f)
		if(exists("datAgOne2.f"))bb=length(datAgOne2.f)
		if(exists("datAgTwo2.f"))cc=length(datAgTwo2.f)
		if(exists("datIllBA2.f"))dd=length(datIllBA2.f)
		if(exists("lumi_NQ.f"))ee=length(lumi_NQ.f)
		if(exists("data.matrix_Nimblegen2.f"))ff=length(data.matrix_Nimblegen2.f)
		if(exists("data.matrixNorm.f"))gg=length(data.matrixNorm.f)
		if(exists("data.matrix_onlineNorm.f"))hh=length(data.matrix_onlineNorm.f)
		}),silent=TRUE)
	if(aa!=0)platforms=c(platforms,"Affymetrix")
	if(bb!=0)platforms=c(platforms,"Agilent_OneColor")
	if(cc!=0)platforms=c(platforms,"Agilent_TwoColor")
	if(dd!=0)platforms=c(platforms,"Illumina_Beadarray")
	if(ee!=0)platforms=c(platforms,"Illumina_Lumi")
	if(ff!=0)platforms=c(platforms,"Nimblegen")
	if(gg!=0)platforms=c(platforms,"Series_Matrix")
	if(hh!=0)platforms=c(platforms,"Online_Data")

	DE_Affy=NULL;DE_Ag1=NULL;DE_Ag2=NULL;DE_Il_B=NULL;DE_Il_L=NULL;DE_N=NULL;DE_S=NULL;DE_O=NULL;
	rm(DE_Affy,DE_Ag1,DE_Ag2,DE_Il_B,DE_Il_L,DE_N,DE_S,DE_O)

	x=NULL
	z=NULL
	w_dge<-gwindow("Select your data",width=260,height=280,visible=FALSE,horizontal=FALSE)
	gp_dge<-ggroup(container=w_dge,horizontal=FALSE)
	cbg_dge<-gcheckboxgroup(platforms,container=gp_dge,handler=f)
	svalue(cbg_dge,index=FALSE)<-1:8

	p_value<-0.01;adjust_value<-"BH";s_value<-"p"
		
	gp_dge_p1<-ggroup(container=gp_dge,horizontal=TRUE)
	numb_l<-gbutton("Number",container=gp_dge_p1,anchor=c(-1,1))
	size(numb_l)=c(80,25)
	numb_numb<-gedit("",initial.msg="Maximum genes to list",width=10,height=20,container=gp_dge_p1,anchor=c(-1,1))
	
	gp_dge_p2<-ggroup(container=gp_dge,horizontal=TRUE)
	logfc_l<-gbutton("logFC",container=gp_dge_p2,anchor=c(-1,1))
	size(logfc_l)=c(80,25)
	logfc_logfc<-gedit("",initial.msg="Log fold change",width=10,height=20,container=gp_dge_p2,anchor=c(-1,1))

	gp_dge_p3<-ggroup(container=gp_dge,horizontal=TRUE)
	pvalue_l<-gbutton("P-value",container=gp_dge_p3,anchor=c(-1,1))
	size(pvalue_l)=c(80,25)
	p_list<-c(0.0001,0.001,0.01,0.05,0.1,0.5,1,10)
	p_value_combo<-gcombobox(p_list,selected=3,container=gp_dge_p3,handler=function(hcp,...){
		x<-svalue(hcp$obj)
		p_value=NULL;
		p_value<<-x
		}
	)
	size(p_value_combo)=c(80,25)
	
	gp_dge_p4<-ggroup(container=gp_dge,horizontal=TRUE)
	adjust_l<-gbutton("Adjustment",container=gp_dge_p4,anchor=c(-1,1))
	size(adjust_l)=c(80,25)
	adjust_list<-c("BH","BY","holm","hochberg","hommel","bonferroni","fdr","none")
	adjust_value_combo<-gcombobox(adjust_list,selected=1,container=gp_dge_p4,handler=function(hcad,...){
		x<-svalue(hcad$obj)
		adjust_value=NULL
		adjust_value<<-x
		}
	)
	size(adjust_value_combo)=c(80,25)
	
	gp_dge_p5<-ggroup(container=gp_dge,horizontal=TRUE)
	sort_l<-gbutton("Sort",container=gp_dge_p5,anchor=c(-1,1))
	size(sort_l)=c(80,25)
	sort_list<-c("B","F","P","p","t","logFC","AveExpr","none")
	p_value_combo<-gcombobox(sort_list,selected=4,container=gp_dge_p5,handler=function(hcso,...){
		x<-svalue(hcso$obj)
		s_value=NULL;
		s_value<<-x
		}
	)
	size(p_value_combo)=c(80,25)

	gp2_dge<-ggroup(container=gp_dge,width=30,height=15,horizontal=TRUE)
	addSpring(gp2_dge)
	y<-gbutton("CANCEL",border=TRUE,handler=function(h,...){
		dispose(w_dge)
		},container=gp2_dge,anchor=c(1,-1))
	y2<-gbutton("OK",border=TRUE,handler=function(h,...){
		dge_n<-svalue(numb_numb)
		if(dge_n=="")dge_n<-10
		dge_lfc<-svalue(logfc_logfc)
		if(dge_lfc=="")dge_lfc=2
		if(length(x)!=0){
			if(length(which(x=="Affymetrix"))!=0){
				err<-try(DE_Affy<<-toptable(dat2Affy.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
				p.value=p_value,sort.by=s_value),silent=TRUE)
				if(dim(DE_Affy)[1]==0 || length(DE_Affy)==0)
				{
					err<-try(DE_Affy<<-toptable(dat2Affy.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Affy)[1]==0 || length(DE_Affy)==0)
				{
					err<-try(DE_Affy<<-toptable(dat2Affy.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Affy)[1]==0 || length(DE_Affy)==0)
				{
					dim_mat<-dim(dat2Affy.f)
					DE_Affy<<-toptable(dat2Affy.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
					sort.by=s_value)
				}
				if(length(DE_Affy)!=0){
					visible(g1_1)<-FALSE
					l$Affymetrix$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
 				}
			if(length(which(x=="Agilent_OneColor"))!=0){
				err<-try(DE_Ag1<<-toptable(datAgOne2.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				sort.by=s_value),silent=TRUE)
				if(dim(DE_Ag1)[1]==0 || length(DE_Ag1)==0)
				{
					err<-try(DE_Ag1<<-toptable(datAgOne2.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Ag1)[1]==0 || length(DE_Ag1)==0)
				{
					err<-try(DE_Ag1<<-toptable(datAgOne2.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Ag1)[1]==0 || length(DE_Ag1)==0)
				{
					dim_mat<-dim(datAgOne2.f)
					DE_Ag1<<-toptable(datAgOne2.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
					sort.by=s_value)
					}
				if(length(DE_Ag1)!=0){
					visible(g1_1)<-FALSE
					l$Agilent_OneColor$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
				}
			if(length(which(x=="Agilent_TwoColor"))!=0){
				err<-try(DE_Ag2<<-toptable(datAgTwo2.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				sort.by=s_value),silent=TRUE)
				if(dim(DE_Ag2)[1]==0 || length(DE_Ag2)==0)
				{
					err<-try(DE_Ag2<<-toptable(datAgTwo2.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Ag2)[1]==0 || length(DE_Ag2)==0)
				{
					err<-try(DE_Ag2<<-toptable(datAgTwo2.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Ag2)[1]==0 || length(DE_Ag2)==0)
				{
					dim_mat<-dim(datAgTwo2.f)
					DE_Ag2<<-toptable(datAgTwo2.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
					sort.by=s_value)
					}
				if(length(DE_Ag2)!=0){
					visible(g1_1)<-FALSE
					l$Agilent_TwoColor$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
				}
			if(length(which(x=="Illumina_Beadarray"))!=0){
				err<-try(DE_Il_B<<-toptable(datIllBA2.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
				p.value=p_value,sort.by=s_value),silent=TRUE)
				if(dim(DE_Il_B)[1]==0 || length(DE_Il_B)==0)
				{
					err<-try(DE_Il_B<<-toptable(datIllBA2.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Il_B)[1]==0 || length(DE_Il_B)==0)
				{
					err<-try(DE_Il_B<<-toptable(datIllBA2.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Il_B)[1]==0 || length(DE_Il_B)==0)
				{
					dim_mat<-dim(datIllBA2.f)
				    DE_Il_B<<-toptable(datIllBA2.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				    sort.by=s_value)
					}
				if(length(DE_Il_B)!=0){
					visible(g1_1)<-FALSE
					l$Illumina_Beadarray$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
				}
			if(length(which(x=="Illumina_Lumi"))!=0){
				err<-try(DE_Il_L<<-toptable(lumi_NQ.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				sort.by=s_value),silent=TRUE)
				if(dim(DE_Il_L)[1]==0 || length(DE_Il_L)==0)
				{
					err<-try(DE_Il_L<<-toptable(lumi_NQ.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Il_L)[1]==0 || length(DE_Il_L)==0)
				{
					err<-try(DE_Il_L<<-toptable(lumi_NQ.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_Il_L)[1]==0 || length(DE_Il_L)==0)
				{
					dim_mat<-dim(lumi_NQ.f)
					DE_Il_L<<-toptable(lumi_NQ.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				sort.by=s_value)
					}
				if(length(DE_Il_L)!=0){
					visible(g1_1)<-FALSE
					l$Illumina_Lumi$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
				}
			if(length(which(x=="Nimblegen"))!=0){
				err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
				p.value=p_value,sort.by=s_value),silent=TRUE)
				if(dim(DE_N)[1]==0 || length(DE_N)==0)
				{
					err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_N)[1]==0 || length(DE_N)==0)
				{
					err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_N)[1]==0 || length(DE_N)==0)
				{
					dim_mat<-dim(data.matrix_Nimblegen2.f)
					DE_N<<-toptable(data.matrix_Nimblegen2.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value)
					}
				if(length(DE_N)!=0){
					visible(g1_1)<-FALSE
					l$Nimblegen$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
				}
			if(length(which(x=="Series_Matrix"))!=0){
				err<-try(DE_S<<-toptable(data.matrixNorm.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
				p.value=p_value,sort.by=s_value),silent=TRUE)
				if(dim(DE_S)[1]==0 || length(DE_S)==0)
				{
					err<-try(DE_S<<-toptable(data.matrixNorm.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_S)[1]==0 || length(DE_S)==0)
				{
					err<-try(DE_S<<-toptable(data.matrixNorm.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_S)[1]==0 || length(DE_S)==0)
				{
					dim_mat<-dim(data.matrixNorm.f)
					DE_S<<-toptable(data.matrixNorm.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				sort.by=s_value)
					}
				if(length(DE_S)!=0){
					visible(g1_1)<-FALSE
					l$Series_Matrix$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
				}
			if(length(which(x=="Online_Data"))!=0){
				err<-try(DE_O<<-toptable(data.matrix_onlineNorm.f,coef=2,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
				p.value=p_value,sort.by=s_value),silent=TRUE)
				if(dim(DE_O)[1]==0 || length(DE_O)==0)
				{
					err<-try(DE_O<<-toptable(data.matrix_onlineNorm.f,coef=3,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_O)[1]==0 || length(DE_O)==0)
				{
					err<-try(DE_O<<-toptable(data.matrix_onlineNorm.f,coef=4,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,
					p.value=p_value,sort.by=s_value),silent=TRUE)
				}
				if(dim(DE_O)[1]==0 || length(DE_O)==0)
				{
					dim_mat<-dim(data.matrix_onlineNorm.f)
					DE_O<<-toptable(data.matrix_onlineNorm.f,number=dge_n,lfc=dge_lfc,adjust.method=adjust_value,p.value=p_value,
				sort.by=s_value)
					}
				if(length(DE_O)!=0){
					visible(g1_1)<-FALSE
					l$Online_Data$DGE<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					}
				display()
  				}
			dispose(w_dge)
			}else{
			gmessage("Plz select the data for Differential Gene Expression","Select Data")
			}
		},container=gp2_dge,anchor=c(1,-1))
	visible(w_dge)<-TRUE
	}
	
