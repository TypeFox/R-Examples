addSkeletonMappings <-
function (mappings,Annotations,points.size.map,min.points.value,max.points.value,points.fill.map,min.points.fill,max.points.fill)
{
	mappings$node.label = list(
		mappingType="passthrough",
		mappingColumn="name",
		mappingColumnType="String",
		visualProperty="NODE_LABEL"
	)
	
	mappings$node.size = list(
		mappingType="continuous",
		mappingColumn=points.size.map,
		mappingColumnType="Double",
		visualProperty="NODE_SIZE",
		points=list(
			point1=list(
				value=max.points.value,
				lesser="95.0",
				equal="95.0",
				greater="95.0"
			),
			point2=list(
				value=min.points.value,
				lesser="35.0",
				equal="35.0",
				greater="35.0"
			)
		)
	)
	
	mappings$node.fill.color = list(
		mappingType="continuous",
		mappingColumn=points.fill.map,
		mappingColumnType="Double",
		visualProperty="NODE_FILL_COLOR",
		points=list(
			point1=list(
				value=min.points.fill,
				lesser="#00FFFF",
				equal="#00FFFF",
				greater="#00FFFF"
			),
			point2=list(
				value=0,
				lesser="#CCCCCC",
				equal="#CCCCCC",
				greater="#CCCCCC"
			),
			point3=list(
				value=max.points.fill,
				lesser="#FF0033",
				equal="#FF0033",
				greater="#FF0033"
			)
		)
	)
	
	if (Annotations)
	{
		mappings$node.border.paint <- list(
			mappingType="discrete",
			mappingColumn="localization",
			mappingColumnType="String",
			visualProperty="NODE_BORDER_PAINT",
			map=list(
				map1=list(
					key="extracellular",
					value="#660066"
				),
				map2=list(
					key="nuclear",
					value="#0000CC"
				),
				map3=list(
					key="cytoplasm",
					value="#990033"
				),
				map4=list(
					key="plasma membrane",
					value="#006600"
				)
			)
		)
	}
	
	return(mappings)
}
