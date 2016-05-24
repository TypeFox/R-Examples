lmgroupedplot <- function(stat.data, map.data, 		# Required -- statistical data; map data
  panel.types, panel.data, 				# Required -- panel types (e.g. 'map', 'labels', 'dot.cl', etc.);
								# 	which columns in dStats to get the data from
  map.link=NULL,						# Specify in a vector (e.g. "c('Subpopulation', 'Poly_Name')") the columns from
								# 	dStats and dMap (respectively) to link the correct map to the correct data line
  nPanels=length(panel.types),				# A count of the number of panels (not meant to actually be specified by a user
								# 	but is used referenced in this function arguments list
  grp.by, cat,						# Required	-- specifies which column in dStats by which to rank and order the output;
								# 	specifies the perceptual groups (e.g. "c(5,5,5,5,2)" to have 4 gorups of 5 followed 
								# 	by a single group of 2			
  colors=brewer.pal(10, "Spectral"),	
  map.color='lightyellow',
  map.all=FALSE,	

  print.file='no', print.res=NA,

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
  plot.grp.spacing=1,
  plot.panel.spacing=1,
  plot.panel.margins=c(0,0,1,0), ...
){					### end function arguments list

mmgroupedplot(stat.data=stat.data, map.data=map.data, 	
  panel.types=panel.types, panel.data=panel.data,
  map.link=map.link,  nPanels=nPanels,
  grp.by=grp.by, cat=cat,  colors=colors,	
  map.color=map.color, map.all=map.all,	
  print.file=print.file, print.res=print.res, panel.att=panel.att,
  plot.header=plot.header, plot.header.size=plot.header.size,  plot.header.color=plot.header.color,
  plot.footer=plot.footer, plot.footer.size=plot.footer.size, plot.footer.color=plot.footer.color,
  plot.width=plot.width, plot.height=plot.height,
  map.spacing=map.spacing, plot.grp.spacing=plot.grp.spacing,
  plot.panel.spacing=plot.panel.spacing, plot.panel.margins=plot.panel.margins
)

}