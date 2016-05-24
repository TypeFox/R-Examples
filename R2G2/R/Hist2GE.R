Hist2GE <-
function(coords, species = 0, grid, goo, nedges, orient, maxAlt = 1e5, colors = "auto", ...){
  ## Prepare input data
  if(species == 0 | length(species) == 0 | length(species) < nrow(coords)) species = rep("nospecies", nrow(coords))

  coords = data.frame(lon = as.numeric(as.character(coords[, 1])), 
		      lat = as.numeric(as.character(coords[, 2])),
		      sp = species)
  coords = coords[is.na(rowSums(coords[,1:2])) == FALSE, ]


  ## retrieve precomputed grid
  precomp = grid

  # convert radian-based grid into lat / lon
  lon = 360 * round((precomp[, 2] / (2 * pi)), 5) - 180 + 0.01
  lat = 180 * round(precomp[, 3] / pi, 5) - 90 + 0.01
  grid = data.frame(lon = lon, lat = lat)


  ## get grid radius and cells spacing
  a.cste = 6378137
  delta = max(diff(precomp[, 3])) * a.cste / 2
  spacing = min(abs(diff(grid[, 1])))


  ## Assign observations to grid cells
  cat("####### Hist2GE\n##### Start assigning observations to cells, please be patient...\n")

  # organise grid and observational data into a common array
  idx.grid = data.frame(numcell = 1:nrow(grid), idx = rep(0, nrow(grid)), lon = grid$lon, lat = grid$lat)
  idx.obs = data.frame(numcell = 1:nrow(coords), idx = rep(1, nrow(coords)), lon = coords[, 1], lat = coords[, 2])
  idx.tot = rbind(idx.grid, idx.obs)

  # sort that array according to lon
  idx.tot = idx.tot[ order(idx.tot[, 3]),]

  # list where our observations are in that array
  targs = which(idx.tot[, 2] == 1)

  # loop to assign observations to grid cells (closest neighbour, OK since grid is consistent)
  ATTR = NULL
  cnt = 0
  pct = 0
  for(i in targs){
    # go at observation i
    cand = idx.tot[i,]

    # slice earth to focus only of subset of cells (i.e. keep the immediate array neighbours of that observation, increases assignment speed)
    left = cand[, 3] - 1.5 * spacing
    right = cand[, 3] + 1.5 * spacing
    upp = cand[, 4] + 1.5 * spacing
    dwn = cand[, 4] + 1.5 * spacing

    slice = idx.tot[ idx.tot[, 3] >= left & idx.tot[, 3] <= right, ]
    
    # find closest cell in grid
    tmp = as.matrix(dist(slice[, 3:4]))
    tmp = tmp[slice[, 1] == cand[, 1] & slice[, 2] == 1 , ]
    tmp = cbind(slice, tmp)
    closest = tmp[tmp[, 2] == 0, ]
    closest = closest[which.min(closest[, 5]), ]

    # assign observation to closest cell
    tmp = c(cand[, 1], closest[, 1])
    ATTR = rbind(ATTR, tmp)

    # print some updates for the user
    if(cnt == (round(length(targs) / 20))){
      pct = pct + 5
      cat("Checked", pct, "% of observations...\n")
      cnt = 0
      }
    cnt = cnt + 1
    }

  cat("##### Compute statistics and produce outputs, please be patient...\n")
  colnames(ATTR) = c("obsnum", "cellnum")
  ATTR = ATTR[order(ATTR[, 1]), ]

  ## Compute species statistics
  final = cbind(coords, ATTR)

  # species occurrences
  occ.cell = table(final[, 5], final[, 3])

  # species richness
  rich.cell = occ.cell
  rich.cell[rich.cell > 0] = 1
  rich.cell = rowSums(rich.cell)

  ## Finalise outputs and 
  grid.focus = grid[rownames(occ.cell),]
  tmp = cbind(as.matrix(grid.focus), as.matrix(rich.cell))
  colnames(tmp) = c("lon", "lat", "NumSpecies")
  out = cbind(tmp, as.matrix(occ.cell))


  ## Producing KML - Species diversity
  kml.name = paste("Hist2GE_Grid_", goo, sep = '')

  Shapes2GE(center = out[,1:2], 
	    nesting = out[,3], 
	    colors = "auto",
	    goo = paste(kml.name, "_divstats.kml", sep = ""),
	    nedges = nedges, 
	    orient = orient, 
	    maxAlt = maxAlt * out[,3], 
	    radius = delta)

  ## Producing KML - Detailled species occurrences
  occ = NULL
  for(i in 4:ncol(out)){
    sp = rep(colnames(out)[i], nrow(out))
    tmp = data.frame(out[, 1:2], sp, out[, i])
    occ = rbind(occ, tmp)
    }
  occ = occ[ occ[, 4] > 0,]

  Shapes2GE(center = occ[, 1:2], 
	    nesting = as.character(occ[, 3]), 
	    colors = "auto",
	    goo = paste(kml.name, "_occurstats.kml", sep = ""),
	    nedges = nedges, 
	    orient = orient, 
	    maxAlt = maxAlt * occ[, 4], 
	    radius = delta)

  # produce output
  cat("##### Done.\n")
  out
  }

