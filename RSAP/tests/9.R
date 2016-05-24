test.exec_infoquery <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    res <- RSAPExecInfoQuery(conn, '0D_NW_M01', '0D_FC_NW_M01_Q0002')
    #str(res)
    #print(res)
    #checkTrue(length(res$D_NW_NETV) >= 10)
    checkTrue(RSAPClose(conn))
}
