Phylo2GE <- function(geo, phy, resol=.1, minAlt = 1e4, maxAlt = 2e6, goo='Phylo2GE.kml'){  

  ## Match phy and geo, make sure that both share their tips and are ordered similarly
  # Note that geo is matched on phy, so geo must contain all the tips of phy, extra taxa of geo are discarded
  geo = geo[match(phy$tip.label, geo[,1]), ]


  ## 1. Get geo and heights for all ancestral nodes
  # 1.1 Compute node geo centroids
  Ntaxa = length(phy$tip.label)
  Ntot = Ntaxa + phy$Nnode

  geo.childs = data.frame(node = 1:nrow(geo), geo[, 2:3])
  
  for(i in Ntot:(Ntaxa + 1)){
    childs = phy$edge[phy$edge[, 1] == i, 2]
    desc = geo.childs[match(childs, geo.childs$node), 2:3]
    if(length(childs) > 2){
      anc.geo = c(i, colMeans(desc))
      } else {
      if(sum(desc[1, ] == desc[2, ]) == 0) { # if the node's descendants have distinct coordinates
	anc.geo = c(i, curvy(.5, desc[1, ], desc[2, ]))
	} else { # if the node's descendants have the same coordinates
	anc.geo = c(i, as.numeric(desc[1, ]))
	}
      }
    geo.childs = rbind(geo.childs, anc.geo)
    }
  gtmp = geo.childs[order(geo.childs[, 1]), 2:3]

  # Compute ancestral node ages (= heights)
  root = phy$edge[1, 1]
  edgenum = nrow(phy$edge)
  ages = NULL
  for (i in 1:edgenum) {
      if (phy$edge[i, 1] == root) 
	  tmp = 0
      else {
	  tmp = ages[which(phy$edge[, 2] == phy$edge[i, 1])]
      }
      ages = c(ages, tmp + phy$edge.length[i])
      }

  tmp = rbind(c(Ntaxa + 1, 0),cbind(phy$edge[, 2], ages))
  tmp = tmp[order(tmp[, 1]), ]


  ## 1.1 Produce meta table, where all infos about nodes are stored.
  # At that stage, we have phy$edge = paths in the tree and meta = positions (XYZ) of nodes
  meta = data.frame(node = tmp[, 1], Lon = gtmp[, 1], Lat = gtmp[, 2], age=tmp[, 2])
  meta$age = 1 - (meta$age / (max(meta$age)))
  meta$age = minAlt + meta$age * maxAlt

  ### 2. Produce Google Earth kml file
  cat("<?xml version=\"1.0\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.0\">\n<Document>",file=goo, append = FALSE)
  cat("<description>Produced using Phylo2GE R script</description>\n<name>R2G2 - Phylo2GE</name>\n<open>0</open>",file=goo, append = TRUE)

  ## 2.1 Loop over all tips and produce corresponding dots
  for(i in 1:Ntaxa){
    xyz = as.numeric(meta[i, 2:4])
    cat("<Placemark>\n<name>", phy$tip.label[i], "</name>\n<LookAt>\n<longitude>", xyz[1], "</longitude>\n<latitude>", xyz[2], "</latitude>\n<range>", xyz[3], "</range>\n</LookAt>", file=goo, append = TRUE)
    cat("<Point>\n<altitudeMode>relativeToGround</altitudeMode>\n<extrude>1</extrude>\n<coordinates>", paste(xyz, collapse=',', sep=' '), "</coordinates>\n</Point>\n</Placemark>\n", file=goo, append = TRUE)
    }
  ## 2.2 Loop over all branches and produce corresponding segments of tree
  cat("<Style id=\"unselectedLine\">\n<LineStyle>\n<color>ff88ffff</color>\n<width>4</width>\n</LineStyle>\n</Style>\n<Style id=\"selectedLine\">\n<LineStyle>\n<color>ff00bbff</color>\n<width>5</width>\n</LineStyle>\n</Style><Folder><name>Edges</name>", file=goo, append = TRUE)

  for(i in 1:nrow(phy$edge)){
    seg = phy$edge[i, ]
    startDD = meta[seg[1], 2:4]
    stopDD = meta[seg[2], 2:4]
    
    # here we draw the vertical bars of the tree (i.e. those uniting edges at a given ancestral node)
    if(sum(startDD[-3] == stopDD[-3]) == 0){ #if the node's descendants have disctinct geographical coords
      tmp = t(sapply(seq(0, 1, by = resol), curvy, startDD, stopDD))
      } else { #if the node's descendants have the same geographical coords
      tmp = matrix(rep(unlist(startDD[1:2]), each = 3), 3, 2)
      }
    pts = cbind(tmp, rep(as.numeric(startDD[3]), nrow(tmp)))
    str = NULL
    for(j in 1:nrow(pts)){
      str = c(str, paste(pts[j, ], collapse = ',', sep = ''))
      }

    final.coo = meta[seg[2], 2:4]
    final.coo = paste(final.coo, collapse = ',')
    cat("<Placemark>\n<styleUrl>#unselectedLine</styleUrl>\n<LineString>\n<altitudeMode>relativeToGround</altitudeMode>\n<coordinates>\n", paste(str, collapse = '\n'), "\n", final.coo, '\n', "</coordinates>\n</LineString>\n</Placemark>\n", file=goo, append = TRUE)
    }
  cat("</Folder></Document>\n</kml>", file=goo, append = TRUE)
  }