require(RnavGraph) || stop("RnavGraph library not available")
require(PairViz) || stop("PairViz library not available")
local({
			data(olive)
			d.olive <- data.frame(olive[,-c(1,2)])
			ng.olive <- ng_data(name = "Olive",
					data = d.olive,
					shortnames = c("p1","p2","s","oleic","l1","l2","a","e"),
					group = as.numeric(olive[,"Area"]),
					labels = as.character(olive[,"Area"])
			)
			
			
			G <- completegraph(shortnames(ng.olive))
			LG <- linegraph(G)
			ng.lg <- ng_graph("3d olive",LG, layout = 'kamadaKawaiSpring')
			ng.lgnot <- ng_graph("4d olive",complement(LG), layout = 'kamadaKawaiSpring')
			
			
			
			nav <- navGraph(ng.olive,
					list(ng.lg,ng.lgnot),
					list(ng_2d(ng.olive,ng.lg,glyphs = eulerian(as(G,"graphNEL"))),ng_2d(ng.olive,ng.lgnot)))
			
			nav1 <- navGraph(ng.olive,
					list(ng.lg,ng.lgnot),
					list(ng_2d(ng.olive,ng.lg,glyphs = eulerian(as(G,"graphNEL"))),ng_2d(ng.olive,ng.lgnot)))
		})
cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_olive.R", package="RnavGraph"),"\n\n\n"))
