addSkeletonDefaults <-
function (defaults,Annotations)
{
	defaults$def.node.color <- list(
		visualProperty = "NODE_FILL_COLOR",
		value = "#CCCC00"
	)
	
	defaults$def.node.selected.paint <- list(
		visualProperty = "NODE_SELECTED_PAINT",
		value = "#CC00CC"
	)
	
	defaults$def.node.shape <- list(
		visualProperty = "NODE_SHAPE", 
		value = "ELLIPSE"
	)
	
	defaults$def.node.size <- list(
		visualProperty = "NODE_SIZE", 
		value = 35.0
	)
	
	defaults$def.node.label.font.face <- list(
		visualProperty = "NODE_LABEL_FONT_FACE", 
		value = "Arial Black,plain,12"
	)
	
	defaults$def.edge.width <- list(
		visualProperty = "EDGE_WIDTH",
		value = 3.0
	)
	
	defaults$def.edge.stroke.selected.paint <- list(
		visualProperty = "EDGE_STROKE_SELECTED_PAINT",
		value = "#999900"
	)
	
	if (Annotations)
	{
		defaults$def.node.border.width <- list(
			visualProperty = "NODE_BORDER_WIDTH",
			value = 7.0
		)
		
		defaults$def.node.border.paint <- list(
			visualProperty = "NODE_BORDER_PAINT",
			value = "SOLID"
		)
		
		defaults$def.node.border.paint <- list(
			visualProperty = "NODE_BORDER_PAINT",
			value = "#999999"
		)
	}
	
	return(defaults)
}
