##' urlTemplate function
##'
##' A function to return a url for a leaflet template for use as map backgrounds with the spplot1 function.
##'
##' Possible templates are: OpenStreetMap_Mapnik, OpenStreetMap_BlackAndWhite, OpenStreetMap_DE, 
##' OpenStreetMap_France, OpenStreetMap_HOT, OpenTopoMap, Thunderforest_OpenCycleMap, Thunderforest_Transport, 
##' Thunderforest_Landscape, Thunderforest_Outdoors, OpenMapSurfer_Roads, OpenMapSurfer_Grayscale, Hydda_Full, 
##' Hydda_Base, MapQuestOpen_OSM, Stamen_Toner, Stamen_TonerBackground, Stamen_TonerLite, Stamen_Watercolor, 
##' Stamen_Terrain, Stamen_TerrainBackground, Stamen_TopOSMRelief, Esri_WorldStreetMap, Esri_WorldTopoMap, 
##' Esri_WorldImagery, Esri_WorldTerrain, Esri_WorldShadedRelief, Esri_WorldPhysical, Esri_OceanBasemap, 
##' Esri_NatGeoWorldMap, Esri_WorldGrayCanvas, Acetate_all, Acetate_terrain, HERE_satelliteDay, 
##' HERE_hybridDayMobile, HERE_hybridDay
##'
##' See \url{http://leaflet-extras.github.io/leaflet-providers/preview/} for other leaflet templates
##'
##' @param name name of the template to use, the default is 'Stamen_Toner'
##' @return url for the leaflet template
##' @export

urlTemplate <- function(name="Stamen_Toner"){

    if(name=="OpenStreetMap_Mapnik"){
        return("http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png")
    }
    else if(name=="OpenStreetMap_BlackAndWhite"){
        return("http://{s}.tiles.wmflabs.org/bw-mapnik/{z}/{x}/{y}.png")
    }
    else if(name=="OpenStreetMap_DE"){
        return("http://{s}.tile.openstreetmap.de/tiles/osmde/{z}/{x}/{y}.png")
    }
    else if(name=="OpenStreetMap_France"){
        return("http://{s}.tile.openstreetmap.fr/osmfr/{z}/{x}/{y}.png")
    }
    else if(name=="OpenStreetMap_HOT"){
        return("http://{s}.tile.openstreetmap.fr/hot/{z}/{x}/{y}.png")
    }
    else if(name=="OpenTopoMap"){
        return("http://{s}.tile.opentopomap.org/{z}/{x}/{y}.png")
    }
    else if(name=="Thunderforest_OpenCycleMap"){
        return("http://{s}.tile.thunderforest.com/cycle/{z}/{x}/{y}.png")
    }
    else if(name=="Thunderforest_Transport"){
        return("http://{s}.tile.thunderforest.com/transport/{z}/{x}/{y}.png")
    }
    else if(name=="Thunderforest_Landscape"){
        return("http://{s}.tile.thunderforest.com/landscape/{z}/{x}/{y}.png")
    }
    else if(name=="Thunderforest_Outdoors"){
        return("http://{s}.tile.thunderforest.com/outdoors/{z}/{x}/{y}.png")
    }
    else if(name=="OpenMapSurfer_Roads"){
        return("http://openmapsurfer.uni-hd.de/tiles/roads/x={x}&y={y}&z={z}")
    }
    else if(name=="OpenMapSurfer_Grayscale"){
        return("http://openmapsurfer.uni-hd.de/tiles/roadsg/x={x}&y={y}&z={z}")
    }
    else if(name=="Hydda_Full"){
        return("http://{s}.tile.openstreetmap.se/hydda/full/{z}/{x}/{y}.png")
    }
    else if(name=="Hydda_Base"){
        return("http://{s}.tile.openstreetmap.se/hydda/base/{z}/{x}/{y}.png")
    }
    else if(name=="MapQuestOpen_OSM"){
        return("http://otile{s}.mqcdn.com/tiles/1.0.0/{type}/{z}/{x}/{y}.{ext}")
    }
    else if(name=="Stamen_Toner"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/toner/{z}/{x}/{y}.png")
    }
    else if(name=="Stamen_TonerBackground"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/toner-background/{z}/{x}/{y}.png")
    }
    else if(name=="Stamen_TonerLite"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/toner-lite/{z}/{x}/{y}.png")
    }
    else if(name=="Stamen_Watercolor"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/watercolor/{z}/{x}/{y}.png")
    }
    else if(name=="Stamen_Terrain"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png")
    }
    else if(name=="Stamen_TerrainBackground"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/terrain-background/{z}/{x}/{y}.png")
    }
    else if(name=="Stamen_TopOSMRelief"){
        return("http://stamen-tiles-{s}.a.ssl.fastly.net/toposm-color-relief/{z}/{x}/{y}.png")
    }
    else if(name=="Esri_WorldStreetMap"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_WorldTopoMap"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_WorldImagery"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_WorldTerrain"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/World_Terrain_Base/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_WorldShadedRelief"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_WorldPhysical"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/World_Physical_Map/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_OceanBasemap"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/Ocean_Basemap/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_NatGeoWorldMap"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Esri_WorldGrayCanvas"){
        return("http://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}")
    }
    else if(name=="Acetate_all"){
        return("http://a{s}.acetate.geoiq.com/tiles/acetate-hillshading/{z}/{x}/{y}.png")
    }
    else if(name=="Acetate_terrain"){
        return("http://a{s}.acetate.geoiq.com/tiles/terrain/{z}/{x}/{y}.png")
    }
    else if(name=="HERE_satelliteDay"){
        return("http://{s}.{base}.maps.cit.api.here.com/maptile/2.1/maptile/{mapID}/satellite.day/{z}/{x}/{y}/256/png8?app_id={app_id}&app_code={app_code}")
    }
    else if(name=="HERE_hybridDayMobile"){
        return("http://{s}.{base}.maps.cit.api.here.com/maptile/2.1/maptile/{mapID}/hybrid.day.mobile/{z}/{x}/{y}/256/png8?app_id={app_id}&app_code={app_code}")
    }
    else if(name=="HERE_hybridDay"){
        return("http://{s}.{base}.maps.cit.api.here.com/maptile/2.1/maptile/{mapID}/hybrid.day/{z}/{x}/{y}/256/png8?app_id={app_id}&app_code={app_code}")
    }
    else{
        stop("Map layer not currently supported. Note that you can ")
    } 

}



# var OpenStreetMap_Mapnik = L.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
#     maxZoom: 19,
#     attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var OpenStreetMap_BlackAndWhite = L.tileLayer('http://{s}.tiles.wmflabs.org/bw-mapnik/{z}/{x}/{y}.png', {
#     maxZoom: 18,
#     attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var OpenStreetMap_DE = L.tileLayer('http://{s}.tile.openstreetmap.de/tiles/osmde/{z}/{x}/{y}.png', {
#     maxZoom: 18,
#     attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var OpenStreetMap_France = L.tileLayer('http://{s}.tile.openstreetmap.fr/osmfr/{z}/{x}/{y}.png', {
#     maxZoom: 19,
#     attribution: '&copy; Openstreetmap France | &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var OpenStreetMap_HOT = L.tileLayer('http://{s}.tile.openstreetmap.fr/hot/{z}/{x}/{y}.png', {
#     maxZoom: 19,
#     attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>, Tiles courtesy of <a href="http://hot.openstreetmap.org/" target="_blank">Humanitarian OpenStreetMap Team</a>'
# });

# var OpenTopoMap = L.tileLayer('http://{s}.tile.opentopomap.org/{z}/{x}/{y}.png', {
#     maxZoom: 16,
#     attribution: 'Map data: &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>, <a href="http://viewfinderpanoramas.org">SRTM</a> | Map style: &copy; <a href="https://opentopomap.org">OpenTopoMap</a> (<a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>)'
# });

# var Thunderforest_OpenCycleMap = L.tileLayer('http://{s}.tile.thunderforest.com/cycle/{z}/{x}/{y}.png', {
#     attribution: '&copy; <a href="http://www.opencyclemap.org">OpenCycleMap</a>, &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var Thunderforest_Transport = L.tileLayer('http://{s}.tile.thunderforest.com/transport/{z}/{x}/{y}.png', {
#     attribution: '&copy; <a href="http://www.opencyclemap.org">OpenCycleMap</a>, &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     maxZoom: 19
# });

# var Thunderforest_Landscape = L.tileLayer('http://{s}.tile.thunderforest.com/landscape/{z}/{x}/{y}.png', {
#     attribution: '&copy; <a href="http://www.opencyclemap.org">OpenCycleMap</a>, &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var Thunderforest_Outdoors = L.tileLayer('http://{s}.tile.thunderforest.com/outdoors/{z}/{x}/{y}.png', {
#     attribution: '&copy; <a href="http://www.opencyclemap.org">OpenCycleMap</a>, &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var OpenMapSurfer_Roads = L.tileLayer('http://openmapsurfer.uni-hd.de/tiles/roads/x={x}&y={y}&z={z}', {
#     maxZoom: 20,
#     attribution: 'Imagery from <a href="http://giscience.uni-hd.de/">GIScience Research Group @ University of Heidelberg</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var OpenMapSurfer_Grayscale = L.tileLayer('http://openmapsurfer.uni-hd.de/tiles/roadsg/x={x}&y={y}&z={z}', {
#     maxZoom: 19,
#     attribution: 'Imagery from <a href="http://giscience.uni-hd.de/">GIScience Research Group @ University of Heidelberg</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var Hydda_Full = L.tileLayer('http://{s}.tile.openstreetmap.se/hydda/full/{z}/{x}/{y}.png', {
#     attribution: 'Tiles courtesy of <a href="http://openstreetmap.se/" target="_blank">OpenStreetMap Sweden</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var Hydda_Base = L.tileLayer('http://{s}.tile.openstreetmap.se/hydda/base/{z}/{x}/{y}.png', {
#     attribution: 'Tiles courtesy of <a href="http://openstreetmap.se/" target="_blank">OpenStreetMap Sweden</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
# });

# var MapQuestOpen_OSM = L.tileLayer('http://otile{s}.mqcdn.com/tiles/1.0.0/{type}/{z}/{x}/{y}.{ext}', {
#     type: 'map',
#     ext: 'jpg',
#     attribution: 'Tiles Courtesy of <a href="http://www.mapquest.com/">MapQuest</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: '1234'
# });

# var MapBox = L.tileLayer('http://{s}.tiles.mapbox.com/v3/MapBox.' + your_api_code + '/{z}/{x}/{y}.png', {
#     attribution: 'Imagery from <a href="http://mapbox.com/about/maps/">MapBox</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd'
# });

# var Stamen_Toner = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/toner/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 0,
#     maxZoom: 20,
#     ext: 'png'
# });

# var Stamen_TonerBackground = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/toner-background/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 0,
#     maxZoom: 20,
#     ext: 'png'
# });

# var Stamen_TonerLite = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/toner-lite/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 0,
#     maxZoom: 20,
#     ext: 'png'
# });

# var Stamen_Watercolor = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/watercolor/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 1,
#     maxZoom: 16,
#     ext: 'png'
# });

# var Stamen_Terrain = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 4,
#     maxZoom: 18,
#     ext: 'png',
#     bounds: [[22, -132], [70, -56]]
# });

# var Stamen_TerrainBackground = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/terrain-background/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 4,
#     maxZoom: 18,
#     ext: 'png',
#     bounds: [[22, -132], [70, -56]]
# });

# var Stamen_TopOSMRelief = L.tileLayer('http://stamen-tiles-{s}.a.ssl.fastly.net/toposm-color-relief/{z}/{x}/{y}.png', {
#     attribution: 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
#     subdomains: 'abcd',
#     minZoom: 0,
#     maxZoom: 20,
#     ext: 'jpg',
#     bounds: [[22, -132], [51, -56]]
# });

# var Esri_WorldStreetMap = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Source: Esri, DeLorme, NAVTEQ, USGS, Intermap, iPC, NRCAN, Esri Japan, METI, Esri China (Hong Kong), Esri (Thailand), TomTom, 2012'
# });

# var Esri_WorldTopoMap = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Esri, DeLorme, NAVTEQ, TomTom, Intermap, iPC, USGS, FAO, NPS, NRCAN, GeoBase, Kadaster NL, Ordnance Survey, Esri Japan, METI, Esri China (Hong Kong), and the GIS User Community'
# });

# var Esri_WorldImagery = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community'
# });

# var Esri_WorldTerrain = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/World_Terrain_Base/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Source: USGS, Esri, TANA, DeLorme, and NPS',
#     maxZoom: 13
# });

# var Esri_WorldShadedRelief = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Source: Esri',
#     maxZoom: 13
# });

# var Esri_WorldPhysical = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/World_Physical_Map/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Source: US National Park Service',
#     maxZoom: 8
# });

# var Esri_OceanBasemap = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/Ocean_Basemap/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic, DeLorme, NAVTEQ, and Esri',
#     maxZoom: 13
# });

# var Esri_NatGeoWorldMap = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; National Geographic, Esri, DeLorme, NAVTEQ, UNEP-WCMC, USGS, NASA, ESA, METI, NRCAN, GEBCO, NOAA, iPC',
#     maxZoom: 16
# });

# var Esri_WorldGrayCanvas = L.tileLayer('http://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}', {
#     attribution: 'Tiles &copy; Esri &mdash; Esri, DeLorme, NAVTEQ',
#     maxZoom: 16
# });



# var Acetate_all = L.tileLayer('http://a{s}.acetate.geoiq.com/tiles/acetate-hillshading/{z}/{x}/{y}.png', {
#     attribution: '&copy;2012 Esri & Stamen, Data from OSM and Natural Earth',
#     subdomains: '0123',
#     minZoom: 2,
#     maxZoom: 18
# });

# var Acetate_terrain = L.tileLayer('http://a{s}.acetate.geoiq.com/tiles/terrain/{z}/{x}/{y}.png', {
#     attribution: '&copy;2012 Esri & Stamen, Data from OSM and Natural Earth',
#     subdomains: '0123',
#     minZoom: 2,
#     maxZoom: 18
# });

# var HERE_satelliteDay = L.tileLayer('http://{s}.{base}.maps.cit.api.here.com/maptile/2.1/maptile/{mapID}/satellite.day/{z}/{x}/{y}/256/png8?app_id={app_id}&app_code={app_code}', {
#     attribution: 'Map &copy; 1987-2014 <a href="http://developer.here.com">HERE</a>',
#     subdomains: '1234',
#     mapID: 'newest',
#     app_id: 'Y8m9dK2brESDPGJPdrvs',
#     app_code: 'dq2MYIvjAotR8tHvY8Q_Dg',
#     base: 'aerial',
#     maxZoom: 20
# });

# var HERE_hybridDayMobile = L.tileLayer('http://{s}.{base}.maps.cit.api.here.com/maptile/2.1/maptile/{mapID}/hybrid.day.mobile/{z}/{x}/{y}/256/png8?app_id={app_id}&app_code={app_code}', {
#     attribution: 'Map &copy; 1987-2014 <a href="http://developer.here.com">HERE</a>',
#     subdomains: '1234',
#     mapID: 'newest',
#     app_id: 'Y8m9dK2brESDPGJPdrvs',
#     app_code: 'dq2MYIvjAotR8tHvY8Q_Dg',
#     base: 'aerial',
#     maxZoom: 20
# });

# var HERE_hybridDay = L.tileLayer('http://{s}.{base}.maps.cit.api.here.com/maptile/2.1/maptile/{mapID}/hybrid.day/{z}/{x}/{y}/256/png8?app_id={app_id}&app_code={app_code}', {
#     attribution: 'Map &copy; 1987-2014 <a href="http://developer.here.com">HERE</a>',
#     subdomains: '1234',
#     mapID: 'newest',
#     app_id: 'Y8m9dK2brESDPGJPdrvs',
#     app_code: 'dq2MYIvjAotR8tHvY8Q_Dg',
#     base: 'aerial',
#     maxZoom: 20
# });