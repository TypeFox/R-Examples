test.read_cube <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    res <- RSAPListCubes(conn)
    print(res)
    #res <- RSAPGetCube(conn, '0D_NW_T01')
    res <- readCube(conn, '0D_NW_T01')
    print(res$INFOOBJECTS)
    res <- RSAPReadCube(conn, '0D_NW_T01', '20120716',chars=list('0D_NW_SORG', '0D_NW_PROD'), kfigures=list('0D_NW_NETV', '0D_NW_QUANT'), options=list(CHANM=list('0D_NW_SORG'),SIGN=list('I'), COMPOP=list('EQ'), LOW=list('1514')))
    str(res)
    print(res)
    checkTrue(length(res$D_NW_NETV) >= 6)
    checkTrue(RSAPClose(conn))
}
