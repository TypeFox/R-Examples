syncSubsample <- function (x, startSearch = min(as.character(x$study.local.timestamp)), 
    endSearch = max(as.character(x$study.local.timestamp)), syncIntervalSecs = 3600, 
    syncAccuracySecs = 60, minEntities = 2, maxEntities = length(unique(x$individual.local.identifier)), 
    mustEntities = NULL, completeSyncsOnly = TRUE, fast = TRUE) 
{
    ## prepare data
    startSearch <- as.numeric(strptime(startSearch, format = "%F %T", 
        tz = "GMT"))
    endSearch <- as.numeric(strptime(endSearch, format = "%F %T", 
        tz = "GMT"))
    time <- as.numeric(strptime(as.character(x$study.local.timestamp), 
        format = "%F %T", tz = "GMT"))
    data <- cbind(x, POSIXct = time, recID = 1:nrow(x))
    
    ## create synchronization events
    syncEvents <- seq(from = startSearch, to = endSearch, by = syncIntervalSecs)
    ## create lower and upper boundary for synchronization events
    lo <- syncEvents - syncAccuracySecs
    hi <- syncEvents + syncAccuracySecs

    ## check for records at syncEvents + - accuracy
    tInLo <- cut(time, c(-Inf, lo, Inf), labels = FALSE, right = FALSE) - 
        1
    tInHi <- cut(time, c(-Inf, hi, Inf), labels = FALSE, right = TRUE)
    syncID <- tInLo
    syncID[tInLo != tInHi] <- NA

    ## get data at sync times
    ## split data by sync times, data not at sync times is omitted
    recsAtSyncEvents <- split(data[, c("individual.local.identifier", 
        "POSIXct", "recID")], syncID)
    ## if no synchronized events present end function
    if (length(recsAtSyncEvents) == 0) {
        return(message("no synchronal events found with this setting for arguments
                       'startSearch', 'endSearch', 'syncIntervalSecs' and 'syncAccuracySecs'"))
    }

    ## only keep one record for each entity at sync time
    ## get created sync IDs
    syncID <- sort(unique(syncID))
    recsAtSyncEvents <- mapply(function(x, y) x[order(abs(y - 
        x$POSIXct)), ], x = recsAtSyncEvents, y = syncEvents[syncID], 
        SIMPLIFY = F)
    recsAtSyncEvents <- lapply(recsAtSyncEvents, function(x) x[!(duplicated(x$individual.local.identifier)), 
        ])

    ## determine all combinations of entities for which to subsample 
    ## get all present combinations of entities at any sync time
    entCombns <- unique(lapply(recsAtSyncEvents, function(x) sort(x$individual.local.identifier)))
    if (!fast) 
        ## create all possible combinations of entities at sync times
        entCombns <- unique(unlist(unlist(lapply(entCombns, function(x1) lapply(min(minEntities, 
            length(x1)):length(x1), function(x2) combn(x1, x2, 
            simplify = FALSE))), recursive = FALSE), recursive = FALSE))
    ## remove combinations with less entities than 'minEntities' or more entities than 'maxEntities'
    entCombns <- entCombns[sapply(entCombns, length) >= minEntities & 
        sapply(entCombns, length) <= maxEntities]
    ## only keep combinations of entities that contain 'mustEntities'
    if (!is.null(mustEntities)) 
        entCombns <- entCombns[sapply(entCombns, function(x) all(mustEntities %in% 
            x))]
    ## if then no combinations of entities left end function
    if (!length(entCombns)) {
        return(message("no synchronal events with this setting for arguments
                       'minEntities', 'maxEntities' and 'mustEntities'"))
    }

    ## for each combination of entities get the respective records at the sync times
    allIn <- lapply(entCombns, function(x1) lapply(recsAtSyncEvents, 
        function(x2) x2$recID[x2$individual.local.identifier %in% 
            x1]))
    ## if 'completeSyncsOnly' is TRUE then remove incomplete synchronization events
    if (completeSyncsOnly) {
        allIn <- lapply(1:length(entCombns), function(x) allIn[[x]][sapply(allIn[[x]], 
            length) >= length(entCombns[[x]])])
    }

    ## create data for overview
    ## get number of entities in each combination
    numberOfEntities <- sapply(entCombns, length)
    ## get number of synchronal events
    numberOfSyncs <- sapply(allIn, length)
    ## get number of pairs of subsequent synchronal events
    numberOfSubsequentSyncs <- sapply(allIn, function(x) sum(diff(as.numeric(names(x))) == 
        1))
    ## get time string of first synchronal event
    firstEvent <- sapply(allIn, function(x) strftime(as.POSIXlt(syncEvents[as.numeric(names(x)[1])], 
        origin = "1970-01-01", tz = "GMT"), format = "%F %T"))
    ## get time string of last synchronal event
    lastEvent <- sapply(allIn, function(x) strftime(as.POSIXlt(syncEvents[as.numeric(names(x)[length(x)])], 
        origin = "1970-01-01", tz = "GMT"), format = "%F %T"))
    ## create overview
    overview <- data.frame(numberOfEntities, numberOfSyncs, numberOfSubsequentSyncs, 
        firstEvent, lastEvent, syncIntervalSecs, syncAccuracySecs, 
        stringsAsFactors = FALSE)

    ## create synchronized subsamples for each combination of entities
    syncData <- lapply(allIn, function(x) data[unlist(x), -((ncol(data) - 
        1):ncol(data))])
    ## get the sync IDs and sync times for the subsamples
    syncID <- lapply(allIn, function(x) rep(as.numeric(names(x)), 
        sapply(x, length)))
    syncTime <- lapply(syncID, function(x) strftime(as.POSIXlt(syncEvents[x], 
        origin = "1970-01-01", tz = "GMT"), format = "%F %T"))
    syncSub <- mapply(function(x1, x2, x3) cbind(x1, syncTime = x2, 
        syncID = x3), x1 = syncData, x2 = syncTime, x3 = syncID, 
        SIMPLIFY = FALSE)
    ## sort each synchronized subsample by sync ID and entity ID
    syncSub <- lapply(syncSub, function(x) x[order(x$syncID, 
        x$individual.local.identifier), ])
    ## sort output
    ordr <- order(overview$numberOfEntities, overview$numberOfSyncs, 
        overview$numberOfSubsequentSyncs, decreasing = TRUE)
    overview <- overview[ordr, ]
    row.names(overview) <- 1:nrow(overview)
    syncSub <- syncSub[ordr]
    entities <- entCombns[ordr]
    
    return(list(overview = overview, data = syncSub, entities = entities))
}
