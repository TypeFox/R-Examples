
odl = list(
    latex=paste(
        'Data by \\href{http://openstreetmap.org}{OpenStreetMap}',
        ' available under the',
        '\\href{http://opendatacommons.org/licenses/odbl}{Open Database License}'
    ),
    markdown=paste(
        'Data by [OpenStreetMap](http://openstreetmap.org) available under the',
        '[Open Database License](http://opendatacommons.org/licenses/odbl)'
    ),
    html=paste(
        'Data by <a href="http://openstreetmap.org">OpenStreetMap</a>,',
        ' available under the',
        '<a href="http://opendatacommons.org/licenses/odbl">Open Database License</a>'
     ),
    text='Data by OpenStreetMap.org available under the Open Database License (opendatacommons.org/licenses/odbl)'
    )

# openstreetmap.org
osm = list(long=list(
      latex=paste(
          ', cartography is licensed as ',
          '\\href{http://creativecommons.org/licenses/by-sa/2.0}{CC BY-SA}.',
          sep=''
      ),
      markdown=paste(
          ', cartography is licensed as [CC BY-SA](http://creativecommons.org/licenses/by-sa/2.0).',
          sep=''
      ),
      html=paste(
          ', cartography is licensed as',
          ' <a href="http://creativecommons.org/licenses/by-sa/2.0">CC BY-SA</a>.',
          sep=''
      ),
      text =paste(
          ', cartography is licensed as CC BY-SA (see www.openstreetmap.org/copyright).',
          sep=''
      )
  ),
  short=list(
      latex='\\copyright \\href{http://openstreetmap.org/copyright}{OpenStreetMap}',
      markdown='&copy; [OpenStreetMap](http://openstreetmap.org/copyright)',
      html= '&copy; <a href="http://openstreetmap.org/copyright">OpenStreetMap</a>',
      text='copyright OpenStreetMap.org'
      )
)



for(D in names(osm$long)){
  osm$long[[D]] = paste(
      osm$short[[D]], 
      ' contributors. ',
      odl[[D]],
      osm$long[[D]], sep=''
  )
}

osmHumanitarian = osm
osmHumanitarian$long$latex = gsub("cartography",
    "cartography by \\\\href{http://hot.openstreetmap.org/about}{Humanitarian OSM team}",
    osmHumanitarian$long$latex)
osmHumanitarian$long$markdown = gsub("cartography",
    "cartography by [Humanitarian OSM team](http://hot.openstreetmap.org/about)",
    osmHumanitarian$long$markdown)
osmHumanitarian$long$text = gsub("cartography",
    "cartography by Humanitarian OSM team (hot.openstreetmap.org)",
    osmHumanitarian$long$text)
osmHumanitarian$long$html = gsub("cartography",
    "cartography by <a href=\"http://hot.openstreetmap.org/about\">Humanitarian OSM team</a>",
    osmHumanitarian$long$html)

osmLandscape = osm
osmLandscape$long$latex = gsub("cartography[[:print:]]+$",
    "cartography by \\\\href{http://www.thunderforest.com/}{Thunderforest}",
    osmLandscape$long$latex)
osmLandscape$long$markdown = gsub("cartography[[:print:]]+$",
    "cartography by [Thunderforest](http://www.thunderforest.com/)",
    osmLandscape$long$markdown)
osmLandscape$long$text = gsub("cartography[[:print:]]+$",
    "cartography by Thunderforest.com",
    osmLandscape$long$text)
osmLandscape$long$html = gsub("cartography[[:print:]]+$",
    "cartography by <a href=\"http://www.thunderforest.com/\">Thunderforest</a>",
    osmLandscape$long$html)

mapquest = mapquestSat = list(
      short=list(
      latex='Tiles courtesy of \\href{http://www.mapquest.com}{MapQuest}',
      text='Tiles courtesy of MapQuest(www.mapquest.com)',
      markdown='Tiles courtesy of [MapQuest](http://www.mapquest.com)',
      html='Tiles courtesy of <a href="http://www.mapquest.com">MapQuest</a>'
    ),
    long=list()
  )
  for(D in names(mapquest$short)){
    mapquest$long[[D]] = paste(
        mapquest$short[[D]],
				odl[[D]],
        sep='. ')
    mapquestSat$long[[D]] = paste(
            mapquest$short[[D]],
        ", portions courtesy NASA/JPL-Caltech and U.S. Depart. of Agriculture, Farm Service Agency.",
        sep='')
  }

  
  stamen = stamenToner = list(
      short=list(
          latex='\\copyright \\href{http://stamen.com}{Stamen Design}',
          markdown='&copy; [OpenStreetMap](http://stamen.com)',
          html= '&copy; <a href="http://stamen.com">Stamen Design</a>',
          text='copyright Stamen Design'
  ),
        long=list(
      html=paste(
          'Map tiles by <a href="http://stamen.com">Stamen Design</a>',
          'under <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>.'
          ),
      latex=paste(
          'Map tiles by \\href{http://stamen.com}{Stamen Design}',
          'under \\href{http://creativecommons.org/licenses/by/3.0}{CC BY 3.0}.'
           ),
      markdown=paste(
          'Map tiles by [Stamen Design](http://stamen.com)',
          'under [CC BY 3.0](http://creativecommons.org/licenses/by/3.0).'
              ),
      text=paste(
          'Map tiles by Stamen Design',
          'under CC BY 3.0 (creativecommons.org/licenses/by/3.0).'
              )
  )
)
for(D in names(stamenToner$long)){
  stamenToner$long[[D]] = paste(
      stamenToner$long[[D]],
      odl[[D]]
      )
}
for(D in names(stamen$long)){
  stamen$long[[D]] = 
      gsub("http://opendatacommons.org/licenses/odbl",
          "http://creativecommons.org/licenses/by/3.0",
          stamenToner$long[[D]])
  stamen$long[[D]] = 
      gsub("Open Database License",
          "CC BY-SA",
          stamen$long[[D]]) 
}

maptoolkit = waze=cartodb = stamenToner
for(D in names(cartodb$long)){
  for(D2 in c('long','short')){
  cartodb[[D2]][[D]] = gsub(
      "Stamen Design", "CartoDB",
          cartodb[[D2]][[D]]
      )
      cartodb[[D2]][[D]] = gsub(
          "stamen.com", "cartodb.com",
              cartodb[[D2]][[D]]
      )
    waze[[D2]][[D]] = gsub(
        "Stamen Design", "Waze mobile",
            waze[[D2]][[D]]
    )
    waze[[D2]][[D]] = gsub(
        "stamen.com", "www.waze.com/legal/notices",
            waze[[D2]][[D]]
    )
    maptoolkit[[D2]][[D]] = gsub(
        "stamen.com", "www.toursprung.com",
        maptoolkit[[D2]][[D]]
    )
    maptoolkit[[D2]][[D]] = gsub(
        "Stamen Design", "Toursprung GmbH",
        maptoolkit[[D2]][[D]]
    )
    
  }
  maptoolkit$long[[D]] = paste(
      maptoolkit$short[[D]],
      odl[[D]]
      )
   waze$long[[D]] = waze$short[[D]]
      
}


esriAttribution = function(name) {
  
  #'http://downloads2.esri.com/ArcGISOnline/docs/tou_summary.pdf'

  long=list(
      'esri'='Esri, HERE, DeLorme, USGS, Intermap, increment P Corp., NRCAN, Esri Japan, METI, Esri China (Hong Kong), Esri (Thailand), TomTom, MapmyIndia, copyright OpenStreetMap contributors, and the GIS User Community',
      'esri-grey' = 'Esri, HERE, DeLorme, MapmyIndia, copyright OpenStreetMap contributors, and the GIS user community',
      'esri-topo' = 'Esri, HERE, DeLorme, TomTom, Intermap, increment P Corp., GEBCO, USGS, FAO, NPS, NRCAN, GeoBase, IGN, Kadaster NL, Ordnance Survey, Esri Japan, METI, Esri China (Hong Kong), swisstopo, MapmyIndia, copyright OpenStreetMap contributors, and the GIS User Community',
      'esri-transport' = 'Esri, HERE, DeLorme, MapmyIndia, copyright OpenStreetMap contributors'
  )
  
  short = 'Esri, OpenStreetMap, and others'
  
  if(!length(grep("^http", name))){
    name = gsub("esri\\.", "esri-", name)
    if(! all(name %in% names(long))){
      warning('name not an esri tile', name)
      return(list())
    }
    name = osmTiles(name)
  }
  
  if(! all(names(name) %in% names(long))){
    warning('name not an esri tile', names(name))
    return(list())
  }
  
  if(!length(grep("esri", names(name)))){
    return(list()) 
  }
  
  weburl = gsub("/tile/?$", '', name)
  



long = long[names(name)]

  list(long=list(
          latex=paste(
              gsub("copyright", "\\\\copyright", long), 
              ' see \\href{',weburl,
              '}{arcgisonline.com}.',
              sep=''
          ),
          markdown=paste(
              gsub("copyright", '&copy;', long), 
              ' see [arcgisonline.com](',
              weburl, ').',
              sep=''
          ),
          html=paste(
              gsub("copyright", '&copy;', long), 
              ' see <a href=\"',
              weburl,
              '\">arcgisonline.com</a>',
              sep=''
          ),
          text =paste(
              long, 
              ' see arcgisonline.com',
              sep=''
          )
      ),
      short=list(
          latex = paste('\\copyright \\href{', weburl,
              '}{', short,'}',sep=''),
          markdown = paste('&copy;','[', short,'](',
              weburl, ')',sep=''),
          html = paste('&copy; <a href=\"', weburl, 
          '\">', short , '</a>',sep=''),
          text=paste('copyright', short)
      )
  )
  
}

openmapAttribution = function(name, type=c('text','latex','markdown','html'), short=FALSE) {

  type = type[1]
  
  if(!is.null(names(name))){
    name = names(name)
  }
  
  shortlong = c('long','short')[short+1]
  
  name = unique(gsub("Red$|Green$|Blue$|Trans$", "", name))
  result = c()
  for(D in name){
        if(length(grep(
                "^osm|opentopomap|openstreetmap|historical|bw.mapnik", 
                D))){        # openstreetmap
          result[D] = osm[[shortlong]][[type]]
        } else if(length(grep("humanitarian",D))){
          result[D] = osmHumanitarian[[shortlong]][[type]]
          
     } else if(length(grep("landscape",D))){
      result[D] = osmLandscape[[shortlong]][[type]]
      
    } else if(length(grep("mapquest|mqcdn",D))){ # mapquest
        if(length(grep("sat/?$",D))){
          result[D] = mapquestSat[[shortlong]][[type]]
        } else {
          result[D] = mapquest[[shortlong]][[type]]
        }
      } else if(length(grep("^waze",D))){ # waze
        result[D] = waze[[shortlong]][[type]]
      } else if(length(grep("maptoolkit",D))){ 
        result[D] = maptoolkit[[shortlong]][[type]]
      } else if(length(grep("stamen",D))){ # mapquest
        if(length(grep("stamen-toner",D))){
          result[D] = stamenToner[[shortlong]][[type]]
        } else {
          result[D] = stamen[[shortlong]][[type]]
        }
    } else if(length(grep("cartodb",D))){ 
      result[D] = cartodb[[shortlong]][[type]]
    } else if(length(grep("^esri",D))){
      result[[D]] =  esriAttribution(D)[[shortlong]][[type]]    
    } else {
      result[D] = NA
    }
} # loop through name
result
}
