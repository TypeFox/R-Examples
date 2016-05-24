get_first_profile <- function(token) {
    id <- suppressWarnings(list_profiles(start.index = 1L, max.results = 1L, token = token)$id)
    if (is.null(id))
        stop("No views (profiles) found on this account.", call. = FALSE)
    return(id)
}

# Get the Anaytics reporting data
#' @importFrom lubridate ymd_h ymd
#' @include query.R
#' @include get-data.R
#' @include profiles.R
get_report <- function(path, query, token, by = NULL) {
    if (is.null(query$profileId)) {
        query$profileId <- get_first_profile(token)
        warning(sprintf("'profileId' was missing. Used first found 'profileId': %s", paste0("ga:", query$profileId)), call. = FALSE)
    }
    if (!grepl("^ga:", query$profileId))
        query$profileId <- paste0("ga:", query$profileId)
    if (is.null(by))
        json_content <- get_data(path, query, token)
    else
        json_content <- fetch_by(path, query, by, token)
    if (is.null(json_content$rows) || length(json_content$rows) == 0) {
        message("No results were obtained.")
        return(invisible(NULL))
    }
    res <- json_content$rows
    # Convert dates to POSIXct with timezone defined in the GA profile
    if (any(grepl("date", names(res), fixed = TRUE))) {
        profile <- json_content$profileInfo
        timezone <- get_profile(profile$accountId, profile$webPropertyId, profile$profileId, token)$timezone
        if (!is.null(res$dateHour))
            res$dateHour <- ymd_h(res$dateHour, tz = timezone)
        if (!is.null(res[["date"]]))
            res[["date"]] <- ymd(res[["date"]], tz = timezone)
        if (!is.null(res$conversionDate))
            res$conversionDate <- ymd(res$conversionDate, tz = timezone)
    }
    attr(res, "profileInfo") <- json_content$profileInfo
    names(json_content$query) <- rename_params(names(json_content$query))
    attr(res, "query") <- json_content$query
    attr(res, "sampled") <- json_content$containsSampledData
    if (!is.null(json_content$containsSampledData) && isTRUE(json_content$containsSampledData)) {
        json_content$sampleSize <- as.numeric(json_content$sampleSize)
        json_content$sampleSpace <- as.numeric(json_content$sampleSpace)
        attr(res, "sampleSize") <- json_content$sampleSize
        attr(res, "sampleSpace") <- json_content$sampleSpace
        samplePerc <- json_content$sampleSize / json_content$sampleSpace * 100
        if (is.null(by))
            warning(sprintf("Data contains sampled data. Used %d sessions (%1.0f%% of sessions). Try to use the 'fetch.by' param to avoid sampling.", json_content$sampleSize, samplePerc), call. = FALSE)
        else
            warning(sprintf("Data contains sampled data. Used %d sessions (%1.0f%% of sessions).", json_content$sampleSize, samplePerc), call. = FALSE)
    }
    return(res)
}
