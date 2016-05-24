test.program_read <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    parms <- list('PROGRAM_NAME' = 'SAPLGRFC')
    res <- RSAPInvoke(conn, "RPY_PROGRAM_READ", parms)
    #str(res$PROG_INF)
    checkEquals('SAPLGRFC', sub("\\s+$", "", res$PROG_INF$PROGNAME))
    checkTrue(length(res$SOURCE_EXTENDED$LINE) > 10)
    checkTrue(RSAPClose(conn))
}
