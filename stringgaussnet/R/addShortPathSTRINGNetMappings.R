addShortPathSTRINGNetMappings <-
function (mappings)
{
	mappings$edge.stroke.unselected.paint = list(
		mappingType="continuous",
		mappingColumn="NIntermediates",
		mappingColumnType="Double",
		visualProperty="EDGE_STROKE_UNSELECTED_PAINT",
		points=list(
			point=list(
				value=0.0,
				lesser="#000000",
				equal="#000000",
				greater="#0000FF"
			)
		)
	)
	
	mappings$edge.line.type = list(
		mappingType="continuous",
		mappingColumn="NIntermediates",
		mappingColumnType="Double",
		visualProperty="EDGE_LINE_TYPE",
		points=list(
			point=list(
				value=0.0,
				lesser="SOLID",
				equal="SOLID",
				greater="LONG_DASH"
			)
		)
	)
	
	mappings$edge.tranparency = list(
		mappingType="continuous",
		mappingColumn="Distance",
		mappingColumnType="Double",
		visualProperty="EDGE_TRANSPARENCY",
		points=list(
			point1=list(
				value=1,
				lesser="500",
				equal="500",
				greater="500"
			),
			point2=list(
				value=10,
				lesser="1",
				equal="1",
				greater="1"
			)
		)
	)
	
	return(mappings)
}
