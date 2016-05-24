# Accumulate probability by grid id ---------------------------------
probPRES<-function(surface, gridvec, detP=1){
  #tmp = as.double(runif(n=(max(gridvec)+1)))
  USE = .C('calc_prob',
                  as.double(surface),                                   #use  =  use surface with prob of at least one wolverine
                  as.integer(gridvec),                                  #grid = vector of grid indices
                  use = as.integer(rep(0,max(gridvec)+1)),              #detection

                  as.double(detP),                                      #detection probability
                  as.integer(length(surface)),                          #pixels
                  as.integer(max(gridvec)),                             #max_grid
                  occP = as.double(runif(n=(max(gridvec)+1)))          #test values, output detection probabilities
                 )$occP[(unique(gridvec)+1)[-1]]
  USE[USE<0]=0
  USE[USE>1]=1
  return(USE)
  }

# Create encounter history -----------------------------------------
encounter.history <- function(map, Parameters, ...){

  add.args <- list(...)
    n_cells <- add.args$n_cells
    grid_layer <- add.args$grid_layer
    filter.map <- add.args$filter.map
    showSteps <- setDefault(add.args$showSteps, F)
    printN <- setDefault(add.args$printN, 0)

  Parameters <- checkParameters(Parameters, add.args)

  n_visits <- Parameters$n_visits
  n_yrs <- Parameters$n_yrs

  lmda <- rep(Parameters$lmda, length.out=(n_yrs - 1))

  if(showSteps == T)
    cat('\nBuilding landscape...\n'); flush.console()

  if(is.null(grid_layer)){
    map <- checkMap(map, filter.map)
    grid_layer <- createGrid(map, Parameters, filter.map)
    n_cells <- length(unique(grid_layer)[unique(grid_layer) > 0])
  }

  encounter_history <- matrix(0, nrow=n_cells, ncol=n_yrs)

  # 1. Place individuals
  if(Parameters$N <= 0 ) stop('Parameters$N <= 0')

  wolv <- addN(Parameters$N, map, Parameters)
  wolv.df <- wolv.dataframe(wolv)

  # 2. Make use surface
  useLayer <- build.useLayer(map, wolv, Parameters)

  # 3. Calculate probability by grid
  P.pres <- matrix(0, nrow=n_cells, ncol=n_yrs)
  P.pres[,1] <- probPRES(useLayer, grid_layer)

  # 3b. (Optional) Output plots with first year data
  if(showSteps){
    printN <- ""  # Cause printN to print to console
    x <- y <- z <- type <- NULL

    prev.ask <- devAskNewPage(ask=T)

    aspRatio <- local({
        if(grepl('longlat', proj4string(map))){
        bbox <- as(extent(map), "SpatialPoints")
        proj4string(bbox) <- CRS(proj4string(map))
        return(unname(mapasp(bbox)))
        } else return(1)
      })

    # Habitat map
    P1 <- ggplot(data=data.frame(coordinates(map)[,1:2], z=getValues(map)),
        aes(x=x, y=y))+geom_raster(aes(fill=z))+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        coord_fixed(ratio=aspRatio)+
        scale_fill_continuous(low='white',high=grey(.7), limits=c(0,1))
    print(P1+labs(fill='Habitat suitability',
                  title='User supplied habitat suitability map'))

    # Grid
    random.colors <- function(n) sample(colors()[!grepl('white|gr.y|black', colors())],n, replace=T)
    print(ggplot(data=data.frame(coordinates(map)[,1:2], z=factor(grid_layer, levels=0:n_cells)),
      aes(x=x, y=y, fill=z))+geom_raster()+
      scale_fill_manual(values=c("white",random.colors(n_cells)),guide=F)+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      coord_fixed(ratio=aspRatio)+
      labs(title='Grid (colored for contrast only)'))

    # Initial locations
    print(P1+geom_point(data=wolv.dataframe(wolv, map), 
      aes(shape=factor(type), colour=factor(type)))+
      scale_colour_hue(l=30)+
        labs(fill='Habitat suitability',
             shape='Type',colour='Type',
             title='Individual activity center locations - Year 1'))


    # Individual use layer example
    ExampleLayer <- build.useLayer(map, list(wolv[[1]][sample(length(wolv[[1]]),1)]), Parameters, Example=T)
    print(ggplot(data.frame(coordinates(map)[,1:2], z=ExampleLayer),
        aes(x=x, y=y, fill=z))+geom_raster()+
      scale_fill_gradientn(colours=c('white',rev(heat.colors(20))))+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      coord_fixed(ratio=aspRatio)+
      labs(fill=' P(use)  \nby pixel',
           title='Probability of use for one individual - example'))

    # Full use surface
    print(ggplot(data.frame(coordinates(map)[,1:2], z=useLayer),
        aes(x=x, y=y, fill=z))+geom_raster()+
      scale_fill_gradientn(colours=c('white',rev(heat.colors(20))))+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      coord_fixed(ratio=aspRatio)+
      labs(fill='Availability\n  by pixel',
           title='Probability of at least one individual present - Year 1'))


    # Gridded use
    print(ggplot(data.frame(coordinates(map)[,1:2],
      z=c(0, P.pres[, 1])[match(grid_layer, unique(c(0,grid_layer)))]),
        aes(x=x, y=y, fill=z))+geom_raster()+
      scale_fill_gradientn(colours=c('white',rev(heat.colors(20))))+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      coord_fixed(ratio=aspRatio)+
      labs(fill='Availability\n   by cell',
           title='Probability of at least one individual present by cell - Year 1'))

    prev.ask <- devAskNewPage(ask=prev.ask)
    rm('prev.ask', 'ExampleLayer')

    cat('\nTotal number of individuals by year\n')
    flush.console()
  }

  # 4. Sample detections in the first year
  encounter_history[, 1] <- sapply(1:n_cells,
    function(i) paste(rbinom(n=n_visits, size=1, prob=P.pres[i,1]), collapse=''))

  # 5. Loop over years to fill in encounter_history
  if(n_yrs > 1){
   for(tt in 2:n_yrs){
      if(printN!=0)
       cat(nrow(wolv.df),' ',sep='',file=printN,append=T)  # Store population sizes by year

      # 6. Calculate population change between t and t+1
      dN <- round(nrow(wolv.df)*(lmda[tt-1] - 1))

      # 7. Implement population change
       if(dN>0){
           new.wolv <- addN(dN, map, Parameters, wolv.df)
           useLayer <- 1-(1-useLayer)*(1-build.useLayer(map, new.wolv, Parameters))
           wolv.df <- rbind(wolv.df, wolv.dataframe(new.wolv))
       } else if(dN<0) {
          lost.wolv <- dropN(abs(dN), map, Parameters, wolv.df)
          useLayer <- 1-(1-useLayer)/(1-build.useLayer(map, lost.wolv[[2]], Parameters))
          useLayer[is.na(useLayer)]<-0
          wolv.df <- wolv.df[-lost.wolv[[1]],]
       }

      # 8. Sample detections and update encounter_history
      P.pres[,tt] <- probPRES(useLayer, grid_layer)

      encounter_history[,tt] <- sapply(1:n_cells,
          function(i) paste(rbinom(n=n_visits, size=1, prob=P.pres[i,tt]), collapse=''))
    }} # End year loop (and if statement)


  if(printN!=0)
    cat(nrow(wolv.df),'\n',sep='',file=printN,append=T)  # Initial population size

  if(showSteps){
    answer <- readline(prompt='\nSave grid to file as raster (y/n)?  ')
    if(tolower(answer)=='y'){
      writeRaster(setValues(map, grid_layer), 'ExampleGrid.tif', overwrite=T)
      cat("  Saved as", paste0(getwd(), '/ExampleGrid.tif'), "\n")
    }

    cat('\nOutputting P.pres and encounter histories by year...\n')
    encounter_history <- data.frame(round(P.pres,3), encounter_history)
      names(encounter_history) <- paste(rep(paste0('Yr', 1:n_yrs),2),
                                      rep(c('p','ch'),each=n_yrs), sep='.')
      encounter_history <- encounter_history[,order(rep(1:n_yrs,2))]
    return(encounter_history)
  } else return(apply(encounter_history,1,paste, collapse=''))
}
