`parconv` <-
function (typ, a2 = NA, c = NA) {
    ifelse(is.na(a2) && !is.na(c), c.a2(typ, c), ifelse(!is.na(a2) && is.na(c), a2.c(typ, a2), NA))
    }

