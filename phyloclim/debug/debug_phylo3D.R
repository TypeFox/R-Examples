loc <- rbind(c(-60.1875, -17.4375),
	c(25.875, 3.9375), 
	c(95.625, 28.6875),
	c(142.875, -20.8125))
	
phy <- structure(list(edge = structure(c(5L, 6L, 6L, 5L, 7L, 7L, 6L, 
1L, 2L, 7L, 3L, 4L), .Dim = c(6L, 2L)), tip.label = c("t2", "t1", 
"t4", "t3"), edge.length = c(0.547188536995525, 0.452811463004475, 
0.452811463004475, 0.278007772293102, 0.721992227706898, 0.721992227706898
), Nnode = 3L), .Names = c("edge", "tip.label", "edge.length", 
"Nnode"), class = "phylo", ploglik = -4.62307622378372, rates = c(0.927342091170134, 
0.99999999, 0.556506257382091, 0.897694723652227, 0.99999999, 
0.77315967803186), message = "relative convergence (4)")

# v.in.ascii in=/Users/Stoffi/grass/phylo3D.txt out=phylo -zn for=standard --o