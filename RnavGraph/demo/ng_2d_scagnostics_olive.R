require(RnavGraph) || stop("RnavGraph library not available")
require(scagnostics) || stop("scagnostics library not available")
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
			
			scagNav(data = ng.olive,
					scags = c("Monotonic", "NotMonotonic",
							"Clumpy", "NotClumpy",
							"Convex", "NotConvex",  
							"Stringy", "NotStringy",
							"Skinny", "NotSkinny",
							"Outlying","NotOutlying",
							"Sparse", "NotSparse",
							"Striated", "NotStriated",
							"Skewed", "NotSkewed"),
					glyphs = hpaths(shortnames(ng.olive)),
					topFrac = 0.15,
					sep = "::")
		})
cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_scagnostics_olive.R", package="RnavGraph"),"\n\n\n"))
