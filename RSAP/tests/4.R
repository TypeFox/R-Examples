test.read_table <- function()
{
    conn <- RSAPConnect("tests/sap.yml")
    parms <- list('DELIMITER' = '|',
                  'ROWCOUNT' = 2,
                  'QUERY_TABLE' = 'T000')
    res <- RSAPInvoke(conn, "RFC_READ_TABLE", parms)
    #str(res$ENTRIES)
    #str(res$DATA)
    checkEquals(2, length(res$DATA$WA))
    
    parms <- list('DELIMITER' = '|',
                  'FIELDS' = list(FIELDNAME = list('CARRID', 'CONNID', 'PRICE', 'SEATSMAX', 'SEATSOCC')),
                  'OPTIONS' = list(TEXT = list("CARRID = 'AA' ", " AND CONNID = 0017 ")),
                  'QUERY_TABLE' = 'SFLIGHTS2')
    res <- RSAPInvoke(conn, "RFC_READ_TABLE", parms)
    #str(res$FIELDS)
    #str(res$DATA)
    checkTrue(length(res$DATA$WA) >= 15)
    checkTrue(RSAPClose(conn))
}
