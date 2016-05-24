lmplot <- function(stat.data, map.data=NULL, 	# Required -- statistical data; map data
  panel.types, panel.data, 			# Required -- panel types (e.g. 'map', 'labels', 'dot.cl', etc.);
							# 	which columns in dStats to get the data from
  map.link=NULL,					# Specify in a vector (e.g. "c('Subpopulation', 'Poly_Name')") the columns from
							# 	dStats and dMap (respectively) to link the correct map to the correct data line
  nPanels=length(panel.types),		# A count of the number of panels (NOT MEANT TO ACTUALLY BE SPECIFIED BY USERS
							# 	but is referenced in this function arguments list
  ord.by, rev.ord=FALSE,						# Column of the stat table by which to rank and order the output
  grouping, 					# Required	-- specifies which column in dStats by which to rank and order the output;
  median.row=FALSE, 					#	specifies the perceptual groups (e.g. "c(5,5,5,5,2)" to have 4 groups of 5 followed 
							# 	by a single group of 2	
  vertical.align='top',

  median.color=gray(.5),	
  colors=brewer.pal(max(grouping), "Spectral"),	

  ## These 3 are legacy arguements that should now 
  #### be specified within the panel att list
  map.all=FALSE, 
  map.color2='lightgray',
  two.ended.maps=FALSE,
  ############
  ############

  print.file='no', print.res=300,
  
  panel.att=vector("list", nPanels),

  ### Options for entire linked micromap plot ###
  plot.header=NA,
  plot.header.size=NA,
  plot.header.color=NA,
  plot.footer=NA,
  plot.footer.size=NA,
  plot.footer.color=NA,
	
  plot.width=7,
  plot.height=7,

  map.spacing=1,
  plot.pGrp.spacing=1,
  plot.panel.spacing=1,
  plot.panel.margins=c(0,0,1,0), ...
){

mmplot(stat.data=stat.data, map.data=map.data, 
  panel.types=panel.types, panel.data=panel.data,
  map.link=map.link, nPanels=nPanels, ord.by=ord.by, rev.ord=rev.ord,
  grouping=grouping, median.row=median.row,
  vertical.align=vertical.align, median.color=median.color,	
  colors=colors,	map.all=map.all, map.color2=map.color2, two.ended.maps=two.ended.maps,
  print.file=print.file, print.res=print.res,
  panel.att=panel.att,
  plot.header=plot.header, plot.header.size=plot.header.size,
  plot.header.color=plot.header.color, plot.footer=plot.footer,
  plot.footer.size=plot.footer.size, plot.footer.color=plot.footer.color,	
  plot.width=plot.width, plot.height=plot.height,
  map.spacing=map.spacing, plot.pGrp.spacing=plot.pGrp.spacing,
  plot.panel.spacing=plot.panel.spacing,
  plot.panel.margins=plot.panel.margins)

}

printLMPlot <- function(plobject, name=NULL, res=300){
 print(plobject, name=name, res=res)
}