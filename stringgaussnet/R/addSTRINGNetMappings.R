addSTRINGNetMappings <-
function (mappings)
{
	mappings$edge.stroke.unselected.paint = list(
		mappingType="discrete",
		mappingColumn="interaction",
		mappingColumnType="String",
		visualProperty="EDGE_STROKE_UNSELECTED_PAINT",
		map=list(
			map1=list(
				key="textmining",
				value="#6666FF"
			),
			map2=list(
				key="cooccurence",
				value="#3333FF"
			),
			map3=list(
				key="knowledge",
				value="#009999"
			),
			map4=list(
				key="coexpression",
				value="#FF0000"
			),
			map5=list(
				key="neighborhood",
				value="#00CCFF"
			),
			map6=list(
				key="experimental",
				value="#009933"
			),
			map7=list(
				key="combined_score",
				value="#000000"
			)
		)
	)
	
	mappings$edge.tranparency = list(
		mappingType="continuous",
		mappingColumn="Score",
		mappingColumnType="Double",
		visualProperty="EDGE_TRANSPARENCY",
		points=list(
			point1=list(
				value=0,
				lesser="1",
				equal="1",
				greater="1"
			),
			point2=list(
				value=1,
				lesser="500",
				equal="500",
				greater="500"
			)
		)
	)
	
	return(mappings)
}
