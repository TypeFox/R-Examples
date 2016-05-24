# \emph{World Federation of Stock Exchanges},
# http://www.world-exchanges.org/statistics


Capitalization <- matrix(c(
    #  2003,     2004,    2005,     2006,     2007,    2008,   Year /
    1328953, 12707578, 3632303, 15421167, 15650832, 9208934, # Euronext US
     888677,  1177517, 1482184,  1700708,  2186550, 1033448, # TSX Group
     585431,   776402,  804014,  1095858,  1298315,  683871, # Australian SE
     278662,   386321,  553073,   818878,  1819100,  647204, # Bombay SE
     714597,   861462, 1054999,  1714953,  2654416, 1328768, # Hong Kong SE
     252893,   363276,  515972,   774115,  1660096,  600281, # NSE India
     360106,   314315,  286190,   917507,  3694348, 1425354, # Shanghai SE
    2953098,  3557674, 4572901,  4614068,  4330921, 3115803, # Tokyo SE
     726243,   940672,  959910,  1322915,  1781132,  948352, # BME Spanish SE
    1079026,  1194516, 1221106,  1637609,  2105197, 1110579, # Deutsche Boerse
    2460064,  2865243, 3058182,  3794310,  3851705, 1868153, # London SE
    2076410,  2441261, 2706803,  3712680,  4222679, 2101745, # Euronext EU
    727103,    826040,  935448,  1212308,  1271047,  857306),# SIX SE
    byrow = TRUE, ncol = 6,
    dimnames = list(
        c("Euronext US", "TSX Group", "Australian SE", "Bombay SE",
            "Hong Kong SE", "NSE India", "Shanghai SE", "Tokyo SE",
            "BME Spanish SE", "Deutsche Boerse", "London SE",
            "Euronext EU", "SIX SE"), 
        as.character(2003:2008)))
