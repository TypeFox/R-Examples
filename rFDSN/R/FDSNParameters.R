### Search parameters list

FDSNParameterList <- function(interface) {
    #Get a list of search parameters, along with their description
    #INPUTS
    #    INTERFACE - Whether you are looking at a station, dataselect, or event query
    #OUTPUTS
    #    PARAMETERS - A list of parameters for the selected interface
    #        $PARAMETER - The parameter used in the query
    #        $TYPE - What kind of data 
    #        $DEFAULT - Default value
    #        $DESCRIPTION - Long description
    #        $EXAMPLE - Example input

    if(!(interface %in% c("station", "dataselect", "event"))) {
        stop("The interface must be either \"station,\" \"dataselect,\" or \"event.\"")
    }
  
    if(interface == "station") {
        params <- c(
            "start", "2001-12-19", "Limit to metadata describing channels operating on or after the specified start time.", "", "date", 
            "end", "2012-12-31", "Limit to metadata describing channels operating on or before the specified end time.", "", "date",
            "startbefore",  2001-12-31,  "Limit to metadata epochs starting before specified time. Applied to channel epochs.", "",  "date",
            "startafter",  2005-12-31,  "Limit to metadata epochs starting after specified time. Applied to channel epochs.", "",  "date",
            "endbefore",  2005-12-31,  "Limit to metadata epochs ending before specified time. Applied to channel epochs.",  "",  "date",
            "endafter",  2005-12-31,  "Limit to metadata epochs ending after specified time. Applied to channel epochs.", "",  "date",
            "net",  "IU",  "Select one or more network codes . Can be SEED codes or data center defined codes.", "",  "string",
            "sta",  "ANMO",  "Select one or more SEED station codes.", "",  "string",
            "loc",  "00",  "Select one or more SEED location identifier. Use -- for \"Blank\" location IDs (ID's containing 2 spaces).",  "",  "string",
            "cha",  "BH1",  "Select one or more SEED channel codes.",  "", "string",
            "minlat",  15.5,  "Southern boundary.",  -90, "float", 
            "maxlat",  25.0,  "Northern boundary.",  90, "float", 
            "minlon",  -170.0,  "Western boundary.",  -180, "float", 
            "maxlon",  170.0,  "Eastern boundary.",  180, "float", 
            "lat",  35,  "Specify the central latitude point.", "", "float", 
            "lon",  170,  "Specify the central longitude point.", "", "float", 
            "maxradius",  20,  "Specify maximum distance from the geographic point defined by latitude and longitude.", "", "float", 
            "minradius",  19,  "Specify minimum distance from the geographic point defined by latitude and longitude.",  0, "float", 
            "level",  "channel",  "Specify level of detail using network, station, channel,or response.", "station",  "string",
            "includerestricted",  TRUE,  "Specify whether results should include information relating to restricted channels.", TRUE, "boolean",
            "includeavailability",  TRUE,  "Specify if results should include information about time series data availability at the channel level.",  FALSE,  "boolean",
            "updatedafter",  "2012-01-01",  "Limit to metadata updated after specified UTC date; updates are data center specific.", "", "date",
            "format",  "text",  "Specify output format. Valid formats include xml and text.",  "xml",  "string",
            "matchtimeseries",  TRUE,  "Retrieve only metadata that is likely to match available time series data",  FALSE,  "boolean",
            "nodata",  404,  "Specify which HTML Status code is returned when no data is found.",  204, "integer")
    }

    if(interface == "dataselect") {
        params <- c(
           "start",  "2010-02-27T06:30:00",  "Specifies the desired start-time for miniSEED data", "", "date",
           "end",  "2010-02-27T10:30:00",  "Specify the end-time for the miniSEED data", "", "date",
           "net", "IU", "Select one or more network codes . Can be SEED codes or data center defined codes.","", "string",
           "sta",  "ANMO",  "Select one or more SEED station codes.",  "",  "string",
           "loc",  "00",  "Select one or more SEED location identifier. Use -- for \"Blank\" location IDs (ID's containing 2 spaces).",  "",  "string",
           "cha",  "BH1",  "Select one or more SEED channel codes.",  "",  "string",
           "quality",  "B",  "Select data based on miniSEED data quality indicator. D, R, Q, M, B. M and B (default) are treated the same and indicate best available.  If M or B are selected, the output data records will be stamped with an M.",  "B", "string", 
           "minimumlength",  0.0,  "Limit results to continuous data segments of a minimum length specified in seconds.",  0.0,  "Float",
           "longestonlly", FALSE,  "Limit results to the longest continuous segment per channel.", FALSE,  "boolean",
           "nodata",  404,  "Specify which HTML Status code is returned when no data is found.",  204, "integer") 
    }

    if(interface == "event") {
        params <- c(
            "start", "2012-11-29", "Limit to events occurring on or after the specified start time.", "", "date",
            "end", "2012-12-01", "Limit to events occurring on or before the specified end time.", "", "date",
            "minlat", 46.8, "Southern boundary.", -90, "float",
            "maxlat", 46.9, "Northern boundary.", 90, "float",
            "minlon", -122, "Western boundary.", -180, "float",
            "maxlon", -121.5, "Eastern boundary.", 180, "float",
            "lat", 40.0, "Specify the central latitude point.", 0.0, "float",
            "lon", 100.0, "Specify the central longitude point.", 0.0, "float",
            "maxradius", 5.0, "Specify maximum distance from the geographic point defined by latitude and longitude.", 180.0, "float",
            "minradius", 1.0, "Specify minimum distance from the geographic point defined by latitude and longitude.", 0.0, "float",
            "mindepth", -1, "Limit to events with depths equal to or greater than the specified depth", "", "float",
            "maxdepth", 20, "Limit to events with depths less than or equal to the specified depth", "", "float",
            "minmag", -1.0, "Limit to events with a magnitude larger than or equal to the specified minimum.", "", "float",
            "maxmag", 8.3, "Limit to events with a magnitude smaller than or equal to the specified maximum.", "", "float",
            "magtype", "M", "Type of Magnitude used to test minimum and maximum limits. Case insensitive. ex. ML Ms mb Mw all preferred", "preferred", "string",
            "catalog", "ISC", "Specify the catalog from which origins and magnitudes will be retrieved (available catalogs)", "NEIC PDE or ISC", "string",
            "contributor", "NEIC PDE-Q", "Limit to events contributed by a specified contributor (available contributors).", "", "string",
            "limit", 25, "Limit the results to the specified number of events", "", "integer",
            "offset", 20, "Return results starting at the event count specified.", 1, "integer",
            "orderby", "time-asc", "Order results by time/ time-asc or magnitude/ magnitude-asc", "", "string",
            "updatedafter", "2000-05-23", "Limit to events updated after the specified time (useful for synchronizing events).", "", "date",
            "includeallorigins", TRUE, "Retrieve all origins or only the primary origin associated with each event.", FALSE, "boolean",
            "includeallmagnitudes", TRUE, "Retrieve all magnitudes for the event, or only the primary magnitude.", FALSE, "boolean",
            "includearrivals", TRUE, "Specify if phase arrivals should be included.", FALSE, "boolean",
            "eventid", 1234, "Retrieve an event based on the unique ID numbers assigned by the IRIS DMC", "", "integer",
            "format", "text", "Specify format. Valid formats include xml and text", "xml", "string",
            "nodata", 404, "Specify which HTML Status code is returned when no data is found.", 204, "integer")
    }

    params.arr <- array(params, dim = c(5, length(params)/5)) 

    parameters <- list(
        parameter = params.arr[1, ],
        type = params.arr[5, ],
        default = params.arr[4, ],
        description = params.arr[3, ],
        example = params.arr[2, ]
    )
} 
