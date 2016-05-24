test.read_table <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    #res <- RSAPReadTable(conn, "SFLIGHTS2")
    res <- readTable(conn, "SFLIGHTS2")
    str(res)
    checkTrue(length(res$CARRID) >= 15)
    res <- RSAPReadTable(conn, "SFLIGHTS2", options=list("CARRID = 'AA' ", " AND CONNID = 0017 "), fields=list('CARRID', 'CONNID', 'PRICE', 'SEATSMAX', 'SEATSOCC'))
    str(res)
    checkTrue(length(res$CARRID) >= 15)
    checkTrue(RSAPClose(conn))
}
