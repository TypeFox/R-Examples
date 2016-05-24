test.changing <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    parms <- list('START_VALUE' = 2,
                  'COUNTER' = 10)
    res <- RSAPInvoke(conn, "STFC_CHANGING", parms)
    #str(res)
    #str(res$RESULT)
    #str(res$COUNTER)
    checkEquals(12, res$RESULT)
    checkEquals(11, res$COUNTER)
    checkTrue(RSAPClose(conn))
}
