test.deep_structure <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    parms <- list('IMPORTSTRUCT' = list(I = 123, C = 'AbCdEf', STR =  'The quick brown fox ...'))
    res <- RSAPInvoke(conn, "STFC_DEEP_STRUCTURE", parms)
    #str(res)
    checkEquals(123, res$ECHOSTRUCT$I)
    checkEquals('AbCdEf', sub("\\s+$", "", res$ECHOSTRUCT$C))
    checkEquals('The quick brown fox ...', res$ECHOSTRUCT$STR)


    parms <- list('IMPORT_TAB' = list(I = list(123, 456), C = list('AbCdEf', 'xyz'), STR =  list('The quick brown fox ...', 'wah wah sauce'), XSTR = list('a', 'b')))
    res <- RSAPInvoke(conn, "STFC_DEEP_TABLE", parms)
    #str(res$EXPORT_TAB)
    checkEquals(123, res$EXPORT_TAB$I[[1]])
    checkEquals(456, res$EXPORT_TAB$I[[2]])
    checkEquals(10, res$EXPORT_TAB$I[[3]])
    checkEquals('AbCdEf', sub("\\s+$", "", res$EXPORT_TAB$C[[1]]))


    checkTrue(RSAPClose(conn))
}
