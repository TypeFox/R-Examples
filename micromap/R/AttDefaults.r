#*** standard attributes used by most panels
#*** these can be altered after being added to the
#*** current panel's attribute list
standard_att <- function(show=FALSE) list(
  panel.header=NA,
  	panel.header.size=1,
	panel.header.color='black',
	panel.header.face='plain',
	panel.header.font=NA,
	panel.header.lineheight=1,


  panel.width=as.numeric(1),
  panel.bgcolor=NA,

  left.margin=NA,
  right.margin=NA,
  panel.margins=c(1, -.25, 1, -1.5),


  graph.grid.major=as.logical(TRUE),
  graph.grid.minor=as.logical(FALSE),
  graph.grid.color='darkgray',
  graph.bgcolor=NA,
  graph.border.color='black',


  xaxis.line.display=as.logical(TRUE),
  xaxis.ticks.display=as.logical(FALSE),
  xaxis.ticks=NA,

  xaxis.labels=NA,
  xaxis.labels.size=1,
  xaxis.labels.angle=NULL,
  xaxis.labels.hjust=NULL,
  xaxis.labels.vjust=NULL,

  xaxis.text.display=as.logical(TRUE),
  xaxis.title=NA,
  

  	xaxis.title.size=1,
	xaxis.title.color='black',
	xaxis.title.face='plain',
	xaxis.title.font=NA,
	xaxis.title.lineheight=1,	
    

  yaxis.line.display=as.logical(FALSE),
  yaxis.line.display=as.logical(FALSE),
  yaxis.ticks.display=as.logical(FALSE),
  yaxis.text.display=as.logical(FALSE),

  yaxis.title=NA,
  yaxis.ticks=NA,
  yaxis.labels=NA
) 




#*** labels ***#
labels_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
			list(text.font=NA,
				text.face='plain',
				text.size=as.numeric(1),
				align='right'))

  tmp.att$xaxis.line.display=as.logical(FALSE)
  tmp.att$xaxis.ticks.display=as.logical(FALSE)
  tmp.att$xaxis.text.display=as.logical(FALSE)
  tmp.att$graph.grid.major=as.logical(FALSE)

  if(show) tmp.att else invisible(tmp.att)
 }



#*** ranks ***#
ranks_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
			list(font=NA,
				face='plain',
				size=as.numeric(1),
				align='right'))

  tmp.att$xaxis.line.display=as.logical(FALSE)
  tmp.att$xaxis.ticks.display=as.logical(FALSE)
  tmp.att$xaxis.text.display=as.logical(FALSE)
  tmp.att$graph.grid.major=as.logical(FALSE)

  if(show) tmp.att else invisible(tmp.att)
 }



#*** dot legend ***#
dot_legend_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
				list(point.size=as.numeric(1.2), 
				  	point.type=as.numeric(19),
					point.border=as.logical(TRUE)))

	tmp.att$xaxis.text.display = FALSE
	tmp.att$xaxis.line.display = FALSE

	tmp.att$graph.border.color = "white"
	tmp.att$graph.grid.major = FALSE


  if(show) tmp.att else invisible(tmp.att)
 }



#*** dot ***#
dot_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
				list(point.size=as.numeric(1), 
				  	point.type=as.numeric(19),
					point.border=as.logical(TRUE),

					median.line=as.logical(FALSE),
					median.line.col='black',
					median.line.typ='longdash',
					median.line.size=1,

					add.line=NA,
					add.line.col='black',
					add.line.typ='longdash',
					add.line.size=1,
					
					connected.dots = F,
	        			connected.col = gray(.6), 
	        			connected.typ = "solid", 
					connected.size = as.numeric(.5)))

  if(show) tmp.att else invisible(tmp.att)
 }


#*** dot.cl ***#
dot_cl_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
				list(point.size=as.numeric(1), 
					point.type=as.numeric(19),
					point.border=as.logical(TRUE),
 
				  	line.width=as.numeric(1),

					median.line=as.logical(FALSE),
					median.line.col='black',
					median.line.typ='longdash',
					median.line.size=1,

					add.line=NA,
					add.line.col='black',
					add.line.typ='longdash',
					add.line.size=1,
					
					connected.dots = F,
	        			connected.col = gray(.6), 
	        			connected.typ = "solid", 
					connected.size = as.numeric(.5)))


  if(show) tmp.att else invisible(tmp.att)
 }


#*** bar ***#
bar_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
				list(graph.bar.size=as.numeric(1),

					median.line=as.logical(FALSE),
					median.line.col='black',
					median.line.typ='longdash',
					median.line.size=1,

					add.line=NA,
					add.line.col='black',
					add.line.typ='longdash',
					add.line.size=1))

  if(show) tmp.att else invisible(tmp.att)
 }


#*** bar.cl ***#
bar_cl_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
				list(graph.bar.size=as.numeric(1),

					median.line=as.logical(FALSE),
					median.line.col='black',
					median.line.typ='longdash',
					median.line.size=1,

					add.line=NA,
					add.line.col='black',
					add.line.typ='longdash',
					add.line.size=1))

  if(show) tmp.att else invisible(tmp.att)
}


#*** box.summary ***#
box_summary_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(), 
				list(graph.bar.size=as.numeric(1),

					median.line=as.logical(FALSE),
					median.line.col='black',
					median.line.typ='longdash',
					median.line.size=1,

					add.line=NA,
					add.line.col='black',
					add.line.typ='longdash',
					add.line.size=1))

  if(show) tmp.att else invisible(tmp.att)
}





#*** map ***#
map_att <- function(show=FALSE) {
  tmp.att <- append(standard_att(),
			list(fill.regions="aggregate",
				map.all=F, 

				outer.hull = F,
				outer.hull.color='black',
				outer.hull.size=1,

				active.border.color='black',
				active.border.size=1,

				inactive.fill='lightgray',
				inactive.border.color=gray(.25),
				inactive.border.size=1,

				withdata.fill='white',
				withdata.border.color=gray(.75),
				withdata.border.size=1,

				median.fill=gray(.5),

				nodata.fill='white',
				nodata.border.color='white',
				nodata.border.size=1
				))

  tmp.att$graph.grid.major=as.logical(FALSE)
  tmp.att$xaxis.line.display=as.logical(FALSE)
  tmp.att$xaxis.ticks.display=as.logical(FALSE)
  tmp.att$xaxis.text.display=as.logical(FALSE)

  if(show) tmp.att else invisible(tmp.att)
 }



sample_att <- function(size=1, type=rep('standard',size), ord.by=NA, grouping=5,
	colors=brewer.pal(max(grouping), "Spectral"), plot.pGrp.spacing=.05,
	plot.panel.margins=c(0,1,0,0), panel.data=list(NA), median.row=FALSE, show=FALSE){

	att <- vector("list", size)
	for(t in 1:size) att[[t]] <- eval(parse(text=paste(type[t],'_att()',sep='')))
	for(t in 1:size) att[[t]] <- append(att[[t]], list(panel.data=NA))
	att <- append(att, list(ord.by=ord.by, grouping=grouping, colors=colors,
			plot.pGrp.spacing=plot.pGrp.spacing, plot.panel.margins=plot.panel.margins, median.row=median.row))	
	return(att)
}



list_att <- function(panel.type) eval(parse(text=paste("print(",panel.type,"_att(show=TRUE))","")))