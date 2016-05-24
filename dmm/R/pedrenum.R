pedrenum <-
function (df) {
# pedrenum()
    df.r <- df
    Id.r <- match(df$Id, df$Id)
    SId.r <- match(df$SId, df$Id)
    DId.r <- match(df$DId, df$Id)
    df.r$Id <- Id.r
    df.r$SId <- SId.r
    df.r$DId <- DId.r
    return(df.r)
}
