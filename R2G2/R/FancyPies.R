FancyPies <-
function(center, obs, nedges = 3, radius = 50000, orient = 0, diag = FALSE){
  ### Make sure that inputs are robust
  if(nedges < 3){
    nedges = 3
    cat("Please use at least nedges = 3\n")
    }
  
  resol = 1
  npts = round(resol / nedges)
  while(npts <= 1){
    resol = resol + 1
    npts = round(resol / nedges)
    }

  ### Preparing data
  # First, we need to locate the area boundaries in the piechart
  #  each area will be marked by an extra point along the edge (herafter area.bound)
  #  hence, each area.bound point will have to computed using curvy and an ad-hoc f values
  #  We first compute these local f values
  
  area.bound = NULL
  if(length(obs) > 1){
    # converting observations to cumulative proportions
    prop = obs / sum(obs)
    prop = c(0, cumsum(prop))

    # assigning proportions to pie edges
    area.edge = as.numeric(cut(prop, nedges))

    # computing edge-specific f values, corresponding to each 
    f.inits = 1 / nedges * (area.edge - 1)
    f.edge = nedges * (prop - f.inits)
    f.edge[length(f.edge)] = 0

    # making a table of this
    area.bound = data.frame(prop, area.edge, f.edge)
    area.bound = area.bound[ 2:(nrow(area.bound) - 1),]
    }

  ### Compute Pie coordinates
  #Get edge limits (startDD and stopDD of ead edge)
  edges.fig = GetEdges(center, radius, nedges, orient)
  edges.fig = edges.fig[nrow(edges.fig):1,]

  #Get detailled edges (using curvy)
  coords = NULL
  for(i in 1:(nrow(edges.fig) - 1)){ #loop over all needed arcs
    # get start / stop points of each arc
    startDD = edges.fig[i, ] 
    stopDD = edges.fig[i + 1, ]
    
    # list of background points to fill each edge
    f.bgrd = seq(0, 1, length.out = npts)
    
    # list of extra points, marking the pie area limits
    edge.nr = i + 1
    f.rel = area.bound[ area.bound[, 2] == i, 3]
    flag = c(rep(0, npts), rep(1, length(f.rel))) #flag these extra points, important for getting area boundaries

    # list of all points
    f.tot = c(f.bgrd, f.rel)

    # computing actual coordinates with curvy
    tmp=t(sapply(f.tot, curvy, startDD, stopDD))

    # storing results in table format, and ordering in increasing position along arc
    tmp = cbind(flag, tmp)
    tmp = tmp[ order(f.tot), ]

    coords = rbind(coords, tmp)
    }
  # little trick to close the pie
  coords[nrow(coords) , 1] = 1

  # number the pie slices (using the flag info)
  sct = cumsum(coords[, 1]) + 1
  coords[, 1] = sct

  # make it clean and neat
  colnames(coords) = c("area", "lonDD", "latDD")

  ### Prepare diag plot (if needed)
  if(diag == TRUE){
    # Identify area limits
    nsct = max(coords[, 1])
    magic = rainbow(nsct - 1)

    plot(coords[, 2:3], type = 'n')
    for(i in 1:(nsct-1)){
      tmp = coords[ coords[, 1] == i , 2:3]
      end = coords[ coords[, 1] == i+1 , 2:3]
      topol = rbind(center, tmp, end, center)
      polygon(topol, col = magic[i])
      }
    }

  # and return output
  coords
  } #end of FancyPies
