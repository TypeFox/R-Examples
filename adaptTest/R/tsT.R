`tsT` <-
function (typ, a = NA, a0 = NA, a1 = NA, a2 = NA) {
    switch (sum(is.na(c(a, a0, a1, a2)))+1,
            NA,
            ifelse(is.na(a), a0a1a2.a(typ,a0,a1,a2), ifelse(is.na(a0), aa1a2.a0(typ,a,a1,a2), ifelse(is.na(a1), aa0a2.a1(typ,a,a0,a2), ifelse(is.na(a2), aa0a1.a2(typ,a,a0,a1), NA)))),
            ifelse(is.na(a1) && is.na(a2), aa0.a1a2(typ,a,a0), NA),
            NA,
            NA
            )
    }

