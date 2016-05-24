#' Check if object is of class bison
#' @param x input
#' @export
is.bison <- function(x) inherits(x, "bison")

#' Check if object is of class bison_solr
#' @param x input
#' @export
is.bison_solr <- function(x) inherits(x, "bison_solr")

bison_blanktheme <- function(){
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = rep(unit(0,"null"),4))
}

# Function to make map form lat/long data
# param x The input
# param geom Geom to use, one of geom_point or geom_jitter.
# param jitter Jitter points.
# param customize Customize ggplot2 plot.
# param ... Further parameters passed on to bisonmap function
bison_map_maker <- function(x, geom, jitter, customize)
{
  long=lat=group=decimalLongitude=decimalLatitude=region=id=NULL
  
  if(is(x, "bison")){ tt <- x$points } else { tt <- x$points }

  # Make lat/long data numeric
  tt$decimalLatitude <- as.numeric(as.character(tt$decimalLatitude))
  tt$decimalLongitude <- as.numeric(as.character(tt$decimalLongitude))
  
  # Remove points that are not physically possible
  tt <- tt[complete.cases(tt$decimalLatitude, tt$decimalLongitude), ]
  tt <- tt[-(which(tt$decimalLatitude <=90 || tt$decimalLongitude <=180)), ]
  
  # Check if points are inside the contintental US, or US+Alaska, or US+Alaska+Hawaii, or even farther
  if(all(all(tt$decimalLatitude < 49) & all(tt$decimalLatitude > 24.7433195) & all(tt$decimalLongitude > -130) & all(tt$decimalLongitude < -66.9513812))){
    # plot only contiguous USA
    contig_layer <- subset(all_states, id != "Alaska" & id != "Hawaii")
    tt_contig <- tt[point.in.polygon(tt$decimalLongitude,tt$decimalLatitude,contig_layer$long,contig_layer$lat)==1,]
    
    p <- ggplot() +
      geom_polygon(aes(x=long, y=lat, group = group), colour = "black", fill = NA, size = 0.25) +
      coord_map(projection="azequalarea") +
      bison_blanktheme()
    p %+% droplevels(subset(all_states, id != "Alaska" & id != "Hawaii")) +
      geom_point(data=droplevels(tt_contig), aes(decimalLongitude, decimalLatitude), size=3, colour = "red", position=jitter) +
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
    
  } else
    if(all(all(tt$decimalLatitude < 75) & all(tt$decimalLatitude > 15) & all(tt$decimalLongitude > -170) & all(tt$decimalLongitude < -66.9513812))){
      # plot contiguous USA + Alaska + Hawaii
      contig_layer <- subset(all_states, id != "Alaska" & id != "Hawaii")
      AK_layer <- subset(all_states, id == "Alaska")
      HI_layer <- subset(all_states, id == "Hawaii")
      tt_contig <- tt[point.in.polygon(tt$decimalLongitude,tt$decimalLatitude,contig_layer$long,contig_layer$lat)==1,]
      tt_AK <- tt[point.in.polygon(tt$decimalLongitude,tt$decimalLatitude,AK_layer$long,AK_layer$lat)==1,]
      tt_HI <- tt[point.in.polygon(tt$decimalLongitude,tt$decimalLatitude,HI_layer$long,HI_layer$lat)==1,]
      
      p <- ggplot() +
        geom_polygon(aes(x=long, y=lat, group = group), colour = "black", fill = NA, size = 0.25) +
        coord_map(projection="azequalarea") +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        bison_blanktheme()
      AK <- p %+% subset(all_states, id == "Alaska") + 
        theme(legend.position = "none") +
        geom_point(data=tt_AK, aes(decimalLongitude, decimalLatitude), size=3, colour = "red", position=jitter)
      HI <- p %+% subset(all_states, id == "Hawaii") + 
        theme(legend.position = "none") +
        geom_point(data=tt_HI, aes(decimalLongitude, decimalLatitude), size=3, colour = "red", position=jitter)
      contiguous <- p %+% subset(all_states, id != "Alaska" & id != "Hawaii") +
        geom_point(data=tt_contig, aes(decimalLongitude, decimalLatitude), size=3, colour = "red", position=jitter)
      thing <- function(aa, bb, cc) {
        grid.newpage()
        vp <- viewport(width = 1.3, height = 1.3, y = 0.6)
        print(aa, vp = vp)
        subvp1 <- viewport(width = 0.4, height = 0.4, x = 0.18, y = 0.25)
        print(bb, vp = subvp1)
        subvp2 <- viewport(width = 0.25, height = 0.25, x = 0.64, y = 0.22)
        print(cc, vp = subvp2)
      }
      thing(contiguous, AK, HI)
      
    } else
    { 
      # plot world (minus Antarticta)
      world <- map_data(map="world") # get world map data
      world <- subset(world, region != "Antarctica") # remove Antarctica
      message("Some of your points are outside the US. Make sure the data is correct")
      ggplot(world, aes(long, lat)) +
        geom_polygon(aes(group=group), fill="white", color="gray40", size=0.2) +
        #           coord_map(projection="mollweide") +
        geom_point(data=tt, aes(decimalLongitude, decimalLatitude), size=3, colour = "red", position=jitter) +
        labs(x="", y="") +
        bison_blanktheme() +
        customize
    }
}

bs_compact <- function (l) Filter(Negate(is.null), l)

mssg <- function(x, y) if(x) message(y)
