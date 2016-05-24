addSIMoNeNetMappings <-
function (mappings)
{
	mappings$edge.stroke.unselected.paint = list(
		mappingType="continuous",
		mappingColumn="Rho",
		mappingColumnType="Double",
		visualProperty="EDGE_STROKE_UNSELECTED_PAINT",
		points=list(
			point1=list(
				value=-1.0,
				lesser="#0000CC",
				equal="#0000CC",
				greater="#0000CC"
			),
			point2=list(
				value=0.0,
				lesser="#999999",
				equal="#999999",
				greater="#999999"
			),
			point3=list(
				value=1.0,
				lesser="#CC0033",
				equal="#CC0033",
				greater="#CC0033"
			)
		)
	)
	
	mappings$edge.width = list(
		mappingType="continuous",
		mappingColumn="P.value",
		mappingColumnType="Double",
		visualProperty="EDGE_WIDTH",
		points=list(
			point1=list(
				value=0.0,
				lesser="20.0",
				equal="20.0",
				greater="20.0"
			),
			point2=list(
				value=0.05,
				lesser="5.0",
				equal="5.0",
				greater="5.0"
			),
			point3=list(
				value=1,
				lesser="1.0",
				equal="1.0",
				greater="1.0"
			)
		)
	)
	
	return(mappings)
}
