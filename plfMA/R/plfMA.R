plfMA <-
function(h,...){
#	library(gWidgets)
#	library(limma)
	options(guiToolkit="RGtk2")
#	require(tcltk)
	try(dispose(wp),silent=TRUE)
	try(y<-h,silent=TRUE)
	if(!exists("y"))h<-matrix(1:100,ncol=10)
	wp<<-gwindow("Quality Control",visible=FALSE,horizontal=TRUE,width=1000,height=600)
	gp1<-gpanedgroup(horizontal=FALSE,spacing=5,use.scrollwindow=FALSE,container=wp)
	gp2<-ggroup(horizontal=TRUE,container=gp1)
	dplot<-gbutton("Scatter-Plot",container=gp2,anchor=c(-1,1));size(dplot)=c(100,40);
	dhist<-gbutton("Histo-gram",container=gp2,anchor=c(-1,1));size(dhist)=c(100,40);
	dbar<-gbutton("Bar-plot",container=gp2,anchor=c(-1,1));size(dbar)=c(100,40);
	dbox<-gbutton("Box-plot",container=gp2,anchor=c(-1,1));size(dbox)=c(100,40);
	dpie<-gbutton("Pie-chart",container=gp2,anchor=c(-1,1));size(dpie)=c(100,40);
	d3d<-gbutton("3D-plot",container=gp2,anchor=c(-1,1));size(d3d)=c(100,40);
	dma<-gbutton("MA-plot",container=gp2,anchor=c(-1,1));size(dma)=c(100,40);
	ddensity<-gbutton("Density-plot",container=gp2,anchor=c(-1,1));size(ddensity)=c(100,40);
	glabel("\t\t",container=gp2)
	import<-gbutton("New",container=gp2,anchor=c(1,-1),spacing=50);size(import)=c(50,25)
	enabled(import)<-TRUE
	export<-gbutton("Export",container=gp2,anchor=c(1,-1));size(export)=c(50,25)
	enabled(export)<-FALSE
	
	gp3<-gpanedgroup(horizontal=TRUE,spacing=20,use.scrollwindow=FALSE,container=gp1,expand=TRUE)
	
	gp4<-ggroup(container=gp3,horizontal=FALSE)
	gp5<-ggroup(container=gp4,horizontal=FALSE)
	gp6<-ggroup(container=gp5,horizontal=FALSE)
	gp7<-ggroup(container=gp6,horizontal=FALSE)

	gp71<-ggroup(container=gp7)
	main<-gbutton("Main",container=gp71,anchor=c(-1,1))
	size(main)=c(50,25)
	main_text<-gedit("",initial.msg="Title",width=27,height=20,container=gp71,anchor=c(-1,1))
	
	gp72<-ggroup(container=gp7)
	sub<-gbutton("Sub",container=gp72,anchor=c(-1,1))
	size(sub)=c(50,25)
	sub_text<-gedit("",initial.msg="Subtitle",width=27,height=20,container=gp72,anchor=c(-1,1))
	
	gp73<-ggroup(container=gp7)
	sub<-gbutton("X-label",container=gp73,anchor=c(-1,1))
	size(sub)=c(50,25)
	x_lab<-gedit("",initial.msg="X-label",width=27,height=20,container=gp73,anchor=c(-1,1))
	
	gp74<-ggroup(container=gp7)
	sub<-gbutton("Y-label",container=gp74,anchor=c(-1,1))
	size(sub)=c(50,25)
	y_lab<-gedit("",initial.msg="Y-label",width=27,height=20,container=gp74,anchor=c(-1,1))

	gp60<-ggroup(container=gp6)
	sub<-gbutton("X-limits",container=gp60,anchor=c(-1,1))
	size(sub)=c(60,25)
	x_lim<-gedit("",initial.msg="xmin,xmax",width=8,height=20,container=gp60,anchor=c(-1,1))
	sub<-gbutton("Y-limits",container=gp60,anchor=c(-1,1))
	size(sub)=c(60,25)
	y_lim<-gedit("",initial.msg="ymin,ymax",width=8,height=20,container=gp60,anchor=c(-1,1))
	
	fc_family=NULL;
	gp600<-ggroup(container=gp6)
	font_family_l<-gbutton("Font Family",container=gp600,anchor=c(-1,1))
	size(font_family_l)=c(100,25)
	families<-c("normal","sans","serif","mono","symbol")
	font_family<-gcombobox(families,selected=1,container=gp600,handler=function(hcf,...){
		x<-svalue(hcf$obj)
		print(x)
		fc_family<<-x
#		fc_family<<-tolower(x)
		print(fc_family)
		}
	)
	size(font_family)=c(150,25)
	
	gp61<-ggroup(container=gp6)
	glabel("\t\t",container=gp61,anchor=c(-1,1))
	col1<-gbutton("Size",container=gp61,anchor=c(-1,1))
	size(col1)=c(68,25)
	col3<-gbutton("Color",container=gp61,anchor=c(-1,1))
	size(col3)=c(68,25)
	col2<-gbutton("Style",container=gp61,anchor=c(-1,1))
	size(col2)=c(68,25)

	f_size<-c("2","3","4","5","6","7","8","9","10","12","14","16","18","20","22","24")
	f_colors<-c("aliceblue","aquamarine","azure","beige","black","blue","brown","chocolate","coral","cornflowerblue",
	"cornsilk","cyan","darkblue","darkcyan","darkgray","darkgreen","darkgrey","darkkhaki","darkmagenta","darkolivegreen",
	"darkorange","darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
	"darkturquoise","darkviolet","deeppink","deepskyblue","dimgray","dimgrey","dodgerblue","firebrick","floralwhite",
	"forestgreen","gainsboro","ghostwhite","gold","goldenrod","gray","green","greenyellow","honeydew","hotpink","indianred",
	"ivory","khaki","lavender","lavenderblush","lawngreen","lemonchiffon","lightblue","lightcoral","lightcyan","lightgoldenrod",
	"lightgoldenrodyellow","lightgray","lightgreen","lightgrey","lightpink","lightsalmon","lightseagreen","lightskyblue",
	"lightslateblue","lightslategray","lightslategrey","lightsteelblue","lightyellow","limegreen","linen","magenta",
	"maroon","mediumaquamarine","mediumblue","mediumorchid","mediumpurple","mediumseagreen","mediumslateblue",
	"mediumspringgreen","mediumturquoise","mediumvioletred","midnightblue","mintcream","mistyrose","moccasin","navajowhite",
	"navy","navyblue","oldlace","olivedrab","orange","orangered","orchid","palegoldenrod","palegreen","paleturquoise",
	"palevioletred","papayawhip","peachpuff","peru","pink","plum","powderblue","purple","red","rosybrown","royalblue",
	"saddlebrown","salmon","sandybrown","seagreen","seashell","sienna","skyblue","slateblue","slategray","slategrey",
	"snow",	"springgreen","steelblue","tan","thistle","tomato","turquoise","violet","violetred","wheat","white","whitesmoke",
	"yellow","yellowgreen")
	f_style<-c("Roman","Bold","Slanted","Bold and Slanted","Symbol")

	fc_sizet=2;fc_colorst="chocolate";fc_stylet=1
	fc_sizes=2;fc_colorss="sienna";fc_styles=1
	fc_sizel=2;fc_colorsl="navy";fc_stylel=1
	fc_sizea=2;fc_colorsa="purple";fc_stylea=1

	gp62<-ggroup(container=gp6)
	titlet<-gbutton("Title",container=gp62,anchor=c(-1,1))
	size(titlet)<-c(65,25)
	font_sizet<-gcombobox(f_size,selected=1,container=gp62,handler=function(hcs,...){
		x<-svalue(hcs$obj)
		print(x)
		fc_sizet<<-x
		}
	)
	size(font_sizet)=c(50,25)
	font_colorst<-gcombobox(f_colors,selected=8,container=gp62,handler=function(hcc,...){
		x<-svalue(hcc$obj)
		print(x)
		fc_colorst<<-x
		}
	)
	size(font_colorst)=c(70,25)
	font_stylet<-gcombobox(f_style,selected=1,container=gp62,handler=function(hcsy,...){
		x<-svalue(hcsy$obj)
		print(x)
		if(x=="Roman")fc_stylet<<-1
		if(x=="Bold")fc_stylet<<-2
		if(x=="Slanted")fc_stylet<<-3
		if(x=="Bold and Slanted")fc_stylet<<-4
		if(x=="Symbol")fc_stylet<<-5
		}
	)
	size(font_stylet)=c(70,25)
	gp63<-ggroup(container=gp6)
	titles<-gbutton("Subtitle",container=gp63,anchor=c(-1,1))
	size(titles)<-c(65,25)
	font_sizes<-gcombobox(f_size,selected=1,container=gp63,handler=function(hcs,...){
		x<-svalue(hcs$obj)
		print(x)
		fc_sizes<<-x
		}
	)
	size(font_sizes)=c(50,25)
	font_colorss<-gcombobox(f_colors,selected=116,container=gp63,handler=function(hcc,...){
		x<-svalue(hcc$obj)
		print(x)
		fc_colorss<<-x
		}
	)
	size(font_colorss)=c(70,25)
	font_styles<-gcombobox(f_style,selected=1,container=gp63,handler=function(hcsy,...){
		x<-svalue(hcsy$obj)
		print(x)
		if(x=="Roman")fc_styles<<-1
		if(x=="Bold")fc_styles<<-2
		if(x=="Slanted")fc_styles<<-3
		if(x=="Bold and Slanted")fc_styles<<-4
		if(x=="Symbol")fc_styles<<-5
		}
	)
	size(font_styles)=c(70,25)

	gp64<-ggroup(container=gp6)
	titlel<-gbutton("Labels",container=gp64,anchor=c(-1,1))
	size(titlel)<-c(65,25)
	font_sizel<-gcombobox(f_size,selected=1,container=gp64,handler=function(hcs,...){
		x<-svalue(hcs$obj)
		print(x)
		fc_sizel<<-x
		}
	)
	size(font_sizel)=c(50,25)
	font_colorsl<-gcombobox(f_colors,selected=90,container=gp64,handler=function(hcc,...){
		x<-svalue(hcc$obj)
		print(x)
		fc_colorsl<<-x
		}
	)
	size(font_colorsl)=c(70,25)
	font_stylel<-gcombobox(f_style,selected=1,container=gp64,handler=function(hcsy,...){
		x<-svalue(hcsy$obj)
		print(x)
		if(x=="Roman")fc_stylel<<-1
		if(x=="Bold")fc_stylel<<-2
		if(x=="Slanted")fc_stylel<<-3
		if(x=="Bold and Slanted")fc_stylel<<-4
		if(x=="Symbol")fc_stylel<<-5
		}
	)
	size(font_stylel)=c(70,25)
	
	gp65<-ggroup(container=gp6)
	titlea<-gbutton("Axis",container=gp65,anchor=c(-1,1))
	size(titlea)<-c(65,25)
	font_sizea<-gcombobox(f_size,selected=1,container=gp65,handler=function(hcs,...){
		x<-svalue(hcs$obj)
		print(x)
		fc_sizea<<-x
		}
	)
	size(font_sizea)=c(50,25)
	font_colorsa<-gcombobox(f_colors,selected=107,container=gp65,handler=function(hcc,...){
		x<-svalue(hcc$obj)
		print(x)
		fc_colorsa<<-x
		}
	)
	size(font_colorsa)=c(70,25)
	font_stylea<-gcombobox(f_style,selected=1,container=gp65,handler=function(hcsy,...){
		x<-svalue(hcsy$obj)
		print(x)
		if(x=="Roman")fc_stylea<<-1
		if(x=="Bold")fc_stylea<<-2
		if(x=="Slanted")fc_stylea<<-3
		if(x=="Bold and Slanted")fc_stylea<<-4
		if(x=="Symbol")fc_stylea<<-5
		}
	)
	size(font_stylea)=c(70,25)
	gp51<-ggroup(container=gp5)
	gp52<-ggroup(container=gp5)
	gp56_1<-ggroup(container=gp5)
	gp56_2<-ggroup(container=gp5)
	gp53<-ggroup(container=gp5)
	gp54<-ggroup(container=gp5)
	gp55<-ggroup(container=gp5)

	gp41<-ggroup(container=gp4)
	glabel("                                          ",container=gp41)
	set<-gbutton("Set",container=gp41,anchor=c(-1,-1))
	size(set)=c(50,25)
	exit<-gbutton("Exit",container=gp41,anchor=c(-1,-1))
	size(exit)=c(50,25)
	
	addHandlerClicked(exit,handler=function(he,...){
		dispose(wp)
		}
	)

	gp8<-ggroup(container=gp3,expand=TRUE)
	plotarea<-ggraphics(ps=5,use.scrollwindow=FALSE,horizontal=FALSE,container=gp8)
	visible(wp)<-TRUE
	bgc<-gbutton("BG color",container=gp51,anchor=c(-1,1))
	size(bgc)=c(63,25)
	fgc<-gbutton("FG color",container=gp51,anchor=c(-1,1))
	size(fgc)=c(63,25)
	addc<-gbutton("Add",container=gp51,anchor=c(-1,1))
	size(addc)=c(63,25)
	bdc<-gbutton("Border",container=gp51,anchor=c(-1,1))
	size(bdc)=c(63,25)

	bgc_colors="cyan";fgc_colors="purple";addc_colors="green";bdc_colors="snow";
	bgc_s<-gcombobox(f_colors,selected=12,container=gp52,handler=function(hbgc,...){
		x<-svalue(hbgc$obj)
		bgc_colors<<-x
		}
	)
	size(bgc_s)<-c(63,25)
	enabled(bgc_s)<-TRUE
	fgc_s<-gcombobox(f_colors,selected=107,container=gp52,handler=function(hfgc,...){
		x<-svalue(hfgc$obj)
		fgc_colors<<-x
		}
	)
	size(fgc_s)<-c(63,25)
	enabled(fgc_s)<-TRUE
	addc_s<-gcombobox(f_colors,selected=44,container=gp52,handler=function(haddc,...){
		x<-svalue(haddc$obj)
		addc_colors<<-x
		}
	)
	size(addc_s)<-c(63,25)
	enabled(addc_s)<-FALSE
	bdc_s<-gcombobox(f_colors,selected=121,container=gp52,handler=function(hbdc,...){
		x<-svalue(hbdc$obj)
		bdc_colors<<-x
		}
	)
	size(bdc_s)<-c(63,25)
	enabled(bdc_s)<-FALSE
	
	pch_1<-gbutton("Point_type",container=gp56_1,anchor=c(-1,1))
	size(pch_1)=c(72,25)
	lty_1<-gbutton("Line_type",container=gp56_1,anchor=c(-1,1))
	size(lty_1)=c(72,25)
	lwd_1<-gbutton("Line_width",container=gp56_1,anchor=c(-1,1))
	size(lwd_1)=c(72,25)
	pch_1_n<-c(1:25,32:127)
	pch_s<-1
	pch_type<-gcombobox(pch_1_n,selected=1,container=gp56_2,handler=function(hpch,...){
		x<-svalue(hpch$obj)
		pch_s<<-x
	}	
	)
	size(pch_type)<-c(72,25)
	lty_1_n<-c("blank","solid","dashed","dotted","dotdash","longdash","twodash")
	lty_s<-"solid"
	lty_type<-gcombobox(lty_1_n,selected=2,container=gp56_2,handler=function(hlty,...){
		x<-svalue(hlty$obj)
		lty_s<<-x
	}
	)
	size(lty_type)<-c(72,25)
	lwd_type<-gedit("",initial.msg="number",width=8,height=20,container=gp56_2,anchor=c(-1,1))
	size(lwd_type)<-c(72,25)
	
	p_ty_l<-gbutton("Type",container=gp53,anchor=c(-1,1))
	size(p_ty_l)=c(65,25)
	types<-c("points","lines","points & lines","lines point apart","lines points overplot","histogram like","stair steps","other steps","none")
	p_ty<-"p"
	plot_type<-gcombobox(types,selected=1,container=gp53,handler=function(htype,...){
		x<-svalue(htype$obj)
		if(x=="points")p_ty<<-"p"
		if(x=="lines")p_ty<<-"l"
		if(x=="points & lines")p_ty<<-"b"
		if(x=="lines point aprt")p_ty<<-"c"
		if(x=="lines points overplot")p_ty<<-"o"
		if(x=="histogram like")p_ty<<-"h"
		if(x=="stair steps")p_ty<<-"s"
		if(x=="other steps")p_ty<<-"S"
		if(x=="none")p_ty<<-"n"
		}
	)
	size(plot_type)=c(150,25)
	enabled(plot_type)<-FALSE

	radius<-gbutton("Radius",container=gp54,anchor=c(-1,1))
	size(radius)=c(65,25)
	types_r<-c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1)
	ty_r<-1
	radius_type<-gcombobox(types_r,selected=1,container=gp54,handler=function(hrtype,...){
		x<-svalue(hrtype$obj)
		ty_r<<-x
		}
	)
	size(radius_type)=c(50,25)
	enabled(radius_type)<-FALSE
	
	shade<-gbutton("Shade",container=gp54,anchor=c(-1,1))
	size(shade)=c(65,25)
	types_sh<-c(0,0.25,0.5,0.75,1)
	ty_sh<-0
	shade_type<-gcombobox(types_sh,selected=1,container=gp54,handler=function(hshtype,...){
		x<-svalue(hshtype$obj)
		ty_sh<<-x
		}
	)
	size(shade_type)=c(50,25)
	enabled(shade_type)<-FALSE

	logb<-gbutton("Log",container=gp55,anchor=c(-1,1))
	size(logb)=c(65,25)
	types_log<-c("","x","y","xy")
	ty_log<-""
	log_type<-gcombobox(types_log,selected=1,container=gp55,handler=function(hlogtype,...){
		x<-svalue(hlogtype$obj)
		ty_log<<-x
		}
	)
	size(log_type)=c(50,25)
	enabled(log_type)<-FALSE
	
	err=NULL;err=plot(h)
	enabled(export)<-TRUE
	
	newdata=NULL
	addHandlerClicked(import,handler=function(h5,...){
		newdata<-tclvalue(tkgetOpenFile())
		if(length(newdata)!=0){
			nd1<-read.table(newdata,sep="\t",header=TRUE)
			h<<-as.matrix(nd1)
			plot(h)
		}
	}
	)

	click1=NULL;click2=NULL;click3=NULL;click4=NULL;click5=NULL;click6=NULL;click7=NULL;click8=NULL;
	addHandlerClicked(dplot,handler=function(h2,...){
		plot(h);
		click1<<-1;click2<<-NULL;click3<<-NULL;click4<<-NULL;click5<<-NULL;click6<<-NULL;click7<<-NULL;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-TRUE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-TRUE
		enabled(bdc_s)<-FALSE
		enabled(addc_s)<-FALSE		
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-TRUE
		}
	)
	addHandlerClicked(dhist,handler=function(h2,...){
		hist(h);
		click1<<-NULL;click2<<-1;click3<<-NULL;click4<<-NULL;click5<<-NULL;click6<<-NULL;click7<<-NULL;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-FALSE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-FALSE
		enabled(bdc_s)<-TRUE
		enabled(addc_s)<-FALSE		
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-FALSE
		}
	)
	
	addHandlerClicked(dbar,handler=function(h2,...){
		barplot(h);
		click1<<-NULL;click2<<-NULL;click3<<-1;click4<<-NULL;click5<<-NULL;click6<<-NULL;click7<<-NULL;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-TRUE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-FALSE
		enabled(bdc_s)<-TRUE
		enabled(addc_s)<-TRUE		
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-TRUE
		}
	)
	
	addHandlerClicked(dbox,handler=function(h2,...){
		boxplot(h);
		click1<<-NULL;click2<<-NULL;click3<<-NULL;click4<<-1;click5<<-NULL;click6<<-NULL;click7<<-NULL;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-TRUE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-FALSE
		enabled(bdc_s)<-TRUE
		enabled(addc_s)<-TRUE		
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-TRUE
		}
	)
	
	addHandlerClicked(dpie,handler=function(h2,...){
		pie(h);
		click1<<-NULL;click2<<-NULL;click3<<-NULL;click4<<-NULL;click5<<-1;click6<<-NULL;click7<<-NULL;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-TRUE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-FALSE
		enabled(bdc_s)<-TRUE
		enabled(addc_s)<-FALSE		
		enabled(radius_type)<-TRUE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-FALSE
		}
	)
	
	addHandlerClicked(d3d,handler=function(h2,...){
		persp(h);
		click1<<-NULL;click2<<-NULL;click3<<-NULL;click4<<-NULL;click5<<-NULL;click6<<-1;click7<<-NULL;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-FALSE
		enabled(y_lim)<-FALSE
		enabled(plot_type)<-FALSE
		enabled(bdc_s)<-TRUE
		enabled(addc_s)<-FALSE		
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-TRUE
		enabled(log_type)<-FALSE
		}
	)
	
	addHandlerClicked(dma,handler=function(h2,...){
		plotMA(h);
		click1<<-NULL;click2<<-NULL;click3<<-NULL;click4<<-NULL;click5<<-NULL;click6<<-NULL;click7<<-1;click8<<-NULL;
		enabled(export)<-TRUE
		enabled(x_lim)<-TRUE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-FALSE
		enabled(bdc_s)<-FALSE
		enabled(addc_s)<-FALSE		
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-TRUE
		}
	)
	
	addHandlerClicked(ddensity,handler=function(h2,...){
		plot(density(h));
		click1<<-NULL;click2<<-NULL;click3<<-NULL;click4<<-NULL;click5<<-NULL;click6<<-NULL;click7<<-NULL;click8<<-1;
		enabled(export)<-TRUE
		enabled(x_lim)<-TRUE
		enabled(y_lim)<-TRUE
		enabled(plot_type)<-TRUE
		enabled(bdc_s)<-FALSE
		enabled(radius_type)<-FALSE
		enabled(shade_type)<-FALSE
		enabled(log_type)<-TRUE
		}
	)


	addHandlerClicked(set,handler=function(h2,...){
		t1<-svalue(main_text)
		t2<-svalue(sub_text)
		t3<-svalue(x_lab)
		t4<-svalue(y_lab)
		t5=NULL;t6=NULL
		if(svalue(x_lim)!=""){
			rx<-svalue(x_lim)
			rx2<-strsplit(rx,",")
			rx3<-c(rx2[[1]][1],rx2[[1]][2])
			t5<-as.numeric(rx3)
			}
		if(svalue(y_lim)!=""){
			ry<-svalue(y_lim)
			ry2<-strsplit(ry,",")
			ry3<-c(ry2[[1]][1],ry2[[1]][2])
			t6<-as.numeric(ry3)
			}
		t7<-as.numeric(svalue(lwd_type))
		if(is.na(t7)==TRUE)t7=1
		par(bg=bgc_colors)
		
#		fc_colorst="chocolate";fc_colorss="sienna";fc_colorsl="navy";fc_colorsa="purple";
#		bgc_colors="cyan";fgc_colors="purple";addc_colors="green";bdc_colors="snow";

		
		if(length(c(click1,click2,click3,click4,click5,click6,click7,click8))==0)
		{
			err=NULL;err=plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Scatter plot such data",title="Plotting Error")
				}
			} else
		
		if(click1==1&&(length(c(click2,click3,click4,click5,click6,click7,click8)))==0)
		{
			err=NULL;err=plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Scatter plot such data",title="Plotting Error")
				}
			} else
		if(click2==1&&(length(c(click1,click3,click4,click5,click6,click7,click8)))==0)
		{
			t5=NULL
			err=NULL;err=hist(h,main=t1,sub=t2,xlab=t3,ylab=t4,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,border=bdc_colors,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Histogram such data",title="Plotting Error")
				}
			} else
		if(click3==1&&(length(c(click1,click2,click4,click5,click6,click7,click8)))==0)
		{
			if(ty_log=="x" || ty_log=="y" || ty_log=="xy")enabled(addc_s)<-FALSE
			if(ty_log=="")enabled(addc_s)<-TRUE
			barplot(h,col.axis=bgc_colors,log=ty_log)
			rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
			err=NULL;err=barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
			add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Barplot such data",title="Plotting Error")
				}
			} else
		if(click4==1&&(length(c(click1,click2,click3,click5,click6,click7,click8)))==0)
		{
			if(ty_log=="x" || ty_log=="y" || ty_log=="xy")enabled(addc_s)<-FALSE
			if(ty_log=="")enabled(addc_s)<-TRUE
			boxplot(h,col.axis=bgc_colors,log=ty_log)
			rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
			err=NULL;err=boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
			add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Boxplot such data",title="Plotting Error")
				}
			} else
		if(click5==1&&(length(c(click1,click2,click3,click4,click6,click7,click8)))==0)
		{		
			err=NULL;err=pie(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,border=bdc_colors,radius=ty_r,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Pie chart such data",title="Plotting Error")
				}
			} else
		if(click6==1&&(length(c(click1,click2,click3,click4,click5,click7,click8)))==0)
		{
			t5=NULL;t6=NULL;
			x=seq(0,1,length.out=nrow(h))
			y=seq(0,1,length.out=ncol(h))
			z=h
			err=NULL;err=persp(x,y,z,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=range(x),ylim=range(y),zlim=range(z,na.rm=TRUE),
			family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,border=bdc_colors,shade=ty_sh,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot 3D plot such data",title="Plotting Error")
				}
			} else
		if(click7==1&&(length(c(click1,click2,click3,click4,click5,click6,click8)))==0)
		{
			err=NULL;err=plotMA(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot MA plot such data",title="Plotting Error")
				}
			} else
		if(click8==1&&(length(c(click1,click2,click3,click4,click5,click6,click7)))==0)
		{
			err=NULL;err=plot(density(h),main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
			cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
			cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
			cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
			cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
			col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
			if(length(grep("Error in",err))!=0){
				gmessage("Cannot Density plot such data",title="Plotting Error")
				}
			}
		}
	)
	addHandlerClicked(export,handler=function(h3,...){
		t1<-svalue(main_text)
		t2<-svalue(sub_text)
		t3<-svalue(x_lab)
		t4<-svalue(y_lab)
		t5=NULL;t6=NULL
		if(svalue(x_lim)!=""){
			rx<-svalue(x_lim)
			rx2<-strsplit(rx,",")
			rx3<-c(rx2[[1]][1],rx2[[1]][2])
			t5<-as.numeric(rx3)
			}
		if(svalue(y_lim)!=""){
			ry<-svalue(y_lim)
			ry2<-strsplit(ry,",")
			ry3<-c(ry2[[1]][1],ry2[[1]][2])
			t6<-as.numeric(ry3)
			}
		t7<-as.numeric(svalue(lwd_type))
		if(is.na(t7)==TRUE)t7=1
		par(bg=bgc_colors)

#		fc_colorst="chocolate";fc_colorss="sienna";fc_colorsl="navy";fc_colorsa="purple";
#		bgc_colors="cyan";fgc_colors="purple";addc_colors="green";bdc_colors="snow";

		if(click1==1&&(length(c(click2,click3,click4,click5,click6,click7,click8)))==0 || length(c(click1,click2,click3,click4,click5,click6,click7,click8))==0)
		{
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						plot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						type=p_ty,col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click2==1&&(length(c(click1,click3,click4,click5,click6,click7,click8)))==0)
		{
			t5=NULL
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						hist(h,main=t1,sub=t2,xlab=t3,ylab=t4,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						hist(h,main=t1,sub=t2,xlab=t3,ylab=t4,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						hist(h,main=t1,sub=t2,xlab=t3,ylab=t4,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						hist(h,main=t1,sub=t2,xlab=t3,ylab=t4,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						hist(h,main=t1,sub=t2,xlab=t3,ylab=t4,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click3==1&&(length(c(click1,click2,click4,click5,click6,click7,click8)))==0)
		{
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=barplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click4==1&&(length(c(click1,click2,click3,click5,click6,click7,click8)))==0)
		{
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=addc_colors)
						err=NULL;err=boxplot(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						add=TRUE,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click5==1&&(length(c(click1,click2,click3,click4,click6,click7,click8)))==0)
		{
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=pie(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,radius=ty_r,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=pie(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,radius=ty_r,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=pie(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,radius=ty_r,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						err=NULL;err=pie(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,radius=ty_r,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						err=NULL;err=pie(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,radius=ty_r,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click6==1&&(length(c(click1,click2,click3,click4,click5,click7,click8)))==0)
		{
			t5=NULL;t6=NULL;
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						x=seq(0,1,length.out=nrow(h))
						y=seq(0,1,length.out=ncol(h))
						z=h
						persp(h)
						err=NULL;err=persp(x,y,z,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=range(x),ylim=range(y),
						zlim=range(z,na.rm=TRUE),family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,shade=ty_sh,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						x=seq(0,1,length.out=nrow(h))
						y=seq(0,1,length.out=ncol(h))
						z=h
						persp(h)
						err=NULL;err=persp(x,y,z,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=range(x),ylim=range(y),
						zlim=range(z,na.rm=TRUE),family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,shade=ty_sh,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						x=seq(0,1,length.out=nrow(h))
						y=seq(0,1,length.out=ncol(h))
						z=h
						persp(h)
						err=NULL;err=persp(x,y,z,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=range(x),ylim=range(y),
						zlim=range(z,na.rm=TRUE),family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,shade=ty_sh,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						x=seq(0,1,length.out=nrow(h))
						y=seq(0,1,length.out=ncol(h))
						z=h
						err=NULL;err=persp(x,y,z,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=range(x),ylim=range(y),
						zlim=range(z,na.rm=TRUE),family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,shade=ty_sh,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						x=seq(0,1,length.out=nrow(h))
						y=seq(0,1,length.out=ncol(h))
						z=h
						err=NULL;err=persp(x,y,z,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=range(x),ylim=range(y),
						zlim=range(z,na.rm=TRUE),family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,shade=ty_sh,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click7==1&&(length(c(click1,click2,click3,click4,click5,click6,click8)))==0)
		{
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=plotMA(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=plotMA(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=plotMA(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						err=NULL;err=plotMA(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						err=NULL;err=plotMA(h,main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,log=ty_log,pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}
		if(click8==1&&(length(c(click1,click2,click3,click4,click5,click6,click7)))==0)
		{
			jpegFileName <- tclvalue(tkgetSaveFile(initialfile = "",
		    filetypes = "{{JPEG} {.jpeg}} {{PNG} {.png}} {{BMP} {.bmp}} {{PDF} {.pdf}} {{POSTSCRIPT} {.postscript}}"
				)
    			)
			if(!nchar(jpegFileName))
			{
			    tkmessageBox(message=paste("The file was not saved",jpegFileName))
			} 
		else{
				
					if(length(grep(".jpeg",jpegFileName))!=0)
					{
						jpeg(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=plot(density(h),main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
 					}
					
				
					else if(length(grep(".png",jpegFileName))!=0)
					{
						png(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=plot(density(h),main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".bmp",jpegFileName))!=0)
					{
						bmp(jpegFileName)
						par(bg=bgc_colors)
						err=NULL;err=plot(density(h),main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
				    }
				    
				
					else if(length(grep(".pdf",jpegFileName))!=0)
					{
						pdf(jpegFileName)
						par(bg=bgc_colors)
						if(length(fc_family)==0 || fc_family=="normal" || fc_family=="symbol")fc_family<-"sans"
						err=NULL;err=plot(density(h),main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				
					else if(length(grep(".postscript",jpegFileName))!=0)
					{
						postscript(jpegFileName)
						par(bg=bgc_colors)
						if(fc_family=="sans" || fc_family=="serif" || fc_family=="mono")fc_family<-"normal"
						err=NULL;err=plot(density(h),main=t1,sub=t2,xlab=t3,ylab=t4,xlim=t5,ylim=t6,family=fc_family,
						cex.main=fc_sizet,col.main=fc_colorst,font.main=fc_stylet,
						cex.sub=fc_sizes,col.sub=fc_colorss,font.sub=fc_styles,
						cex.lab=fc_sizel,col.lab=fc_colorsl,font.lab=fc_stylel,
						cex.axis=fc_sizea,col.axis=fc_colorsa,font.axis=fc_stylea,
						col=fgc_colors,bg=bgc_colors,border=bdc_colors,log=ty_log,
						pch=pch_s,lty=lty_s,lwd=t7)
						dev.off()
    				}
				    
				tkmessageBox(message=paste("The file was saved",jpegFileName))
				}
			}

		}
	)
}
