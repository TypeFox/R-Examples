test.get_table <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    parms <- list('BYPASS_BUFFER' = 'X',
                  'MAX_ENTRIES' = 10,
                  'TABLE_NAME' = 'T005')
    res <- RSAPInvoke(conn, "RFC_GET_TABLE_ENTRIES", parms)
    #str(res$ENTRIES)
    checkEquals(10, length(res$ENTRIES$WA))
    checkTrue(RSAPClose(conn))
}
