`vol.m3` <-
function(dbh.cm, height.m, multiplier=0.5) {
    vol.m3 <- pi * (dbh.cm/200)^2 * height.m * multiplier
    }

