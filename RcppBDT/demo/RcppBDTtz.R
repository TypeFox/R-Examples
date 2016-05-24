
demo.RcppBDTtz  <- function() {

    require(utils, quiet=TRUE, warn=FALSE)
    require(Rcpp, quiet=TRUE, warn=FALSE)
    require(RcppBDT, quiet=TRUE, warn=FALSE)

    tz <- new(bdtTz, "America/Chicago")

    # format() as cat() enforces txt
    cat("tz object initialized as:       ", format(tz), "\n")
    cat("zone std and dst abbreviation:  ", tz$getStdZoneAbbrev(), "and", tz$getDstZoneAbbrev(), "\n")
    cat("zone std and dst names:         ", tz$getStdZoneName(), "and", tz$getDstZoneName(), "\n")

    cat("\n")
    # format() as cat() enforces txt
    cat("2012 year dst start:            ", format(tz$getDstLocalStart(2012)), "\n")
    cat("2012 year dst end:              ", format(tz$getDstLocalEnd(2012)), "\n")

    cat("\n")
    cat("Offset to UTC (in seconds):     ", tz$getUtcOffset(), "\n")
    cat("DST offset (in seconds):        ", tz$getDstOffset(), "\n")

    cat("\n")
    cat("Formal POSIX string of region:  ", tz$getPosixString(), "\n")

    cat("\n")
    tz1 <- new(bdtTz, "America/Bogota")
    tz2 <- new(bdtTz, "America/Belize")
    cat("UTC difference between Bogota and Belize:    ", tz1$getUtcOffset() - tz2$getUtcOffset(), "\n")

    tz1 <- new(bdtTz, "Australia/Sydney")
    tz2 <- new(bdtTz, "Asia/Singapore")
    cat("UTC difference between Sidney and Singapore: ", tz1$getUtcOffset() - tz2$getUtcOffset(), "\n")

    tz1 <- new(bdtTz, "Europe/Moscow")
    tz2 <- new(bdtTz, "Europe/Madrid")
    cat("UTC difference between Moscow and Madrid:    ", tz1$getUtcOffset() - tz2$getUtcOffset(), "\n")

    cat("\n")
    cat("First 75 tz regions in database:\n")
    print(head(tz$getAllRegions(), 75))


}

demo.RcppBDTtz()
