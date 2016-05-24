test.connecting <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    checkTrue(typeof(attr(conn, 'handle_ptr')) == 'externalptr')
    checkTrue(!is.null(attr(conn, 'handle_ptr')))
    info <- RSAPGetInfo(conn)
    #str(info)
    checkEquals(info[['partnerHost']], "nplhost")
    checkEquals(info[['sysNumber']], "42")
    checkEquals(info[['language']], "E")
    checkEquals(info[['progName']], "SAPLSRFC")
    checkEqualsNumeric(as.numeric(info[['sysNumber']]), 42)
    checkTrue(RSAPClose(conn))

    conn <- RSAPConnect(ashost="nplhost", sysnr="42",
                        client="001", user="developer", 
                        passwd="developer", lang="EN", 
                        trace="1", lcheck="1")
    checkTrue(typeof(attr(conn, 'handle_ptr')) == 'externalptr')
    checkTrue(!is.null(attr(conn, 'handle_ptr')))
    info <- RSAPGetInfo(conn)
    #str(info)
    checkEquals(info[['partnerHost']], "nplhost")
    checkEquals(info[['sysNumber']], "42")
    checkEquals(info[['language']], "E")
    checkEquals(info[['progName']], "SAPLSRFC")
    checkEquals(info[['rfcRole']], "C")
    checkEqualsNumeric(as.numeric(info[['sysNumber']]), 42)
    checkTrue(RSAPClose(conn))
}
           
#test.deactivation <- function()
#{
# DEACTIVATED('Deactivating this test function')
#}
