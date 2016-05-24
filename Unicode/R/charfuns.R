u_char_decomposition <-
function(x)
{
    ## For now, this really gives the Decomposition_Mapping property.
    ## Ideally, this should have a type argument controlling whether to
    ## perform canonical or compatibility decomposition, and most likely
    ## also an argument controlling whether to decompose recursively or
    ## not.

    ## Maybe return a u_char_seq eventually?
    
    x <- as.u_char(x)
    ## Could make this more efficient eventually ...
    p <- match(x, UCD_Unicode_data_table$Code)
    y <- sub(".*> *", "", UCD_Unicode_data_table$Decomposition[p])
    p <- which(is.na(p))
    if(length(p)) {
        r <- u_char_match(x[p], UCD_Unicode_data_range$Range, 0L)
        ind <- (r == UCD_Unicode_data_range_pos_Hangul)
        if(any(ind)) {
            pos <- p[ind]
            y[pos] <-
                gsub("U+", "",
                     sapply(unclass(.u_char_decomposition_Hangul(x[pos])),
                            paste, collapse = " "),
                     fixed = TRUE)
        }
        ind <- (r > 0L) & (r != UCD_Unicode_data_range_pos_Hangul)
        if(any(ind)) {
            y[p[ind]] <- ""
        }
    }
    y
}

.u_char_decomposition_Hangul <-
function(s)
{
    y <- rep.int(list(integer()), length(s))
    s <- s - Jamo_S_base
    ok <- (s >= 0) & (s < Jamo_S_count)
    if(any(ok)) {
        s <- s[ok]
        L <- Jamo_L_base + s %/% Jamo_N_count
        V <- Jamo_V_base + (s %% Jamo_N_count) %/% Jamo_T_count
        T <- Jamo_T_base + s %% Jamo_T_count
        y[ok] <- Map("c", L, V, ifelse(T == Jamo_T_base, list(NULL), T))
    }
    as.u_char_seq(y)
}

u_char_from_name <-
function(x, type = c("exact", "grep"), ...)
{
    x <- as.character(x)
    type <- match.arg(type)
    if(type == "exact") {
        ## Matching could/should be less strict ...    
        p <- match(x, UCD_Unicode_data_table$Name)
        y <- UCD_Unicode_data_table$Code[p]
        if(any(is.na(p))) {
            re <- "^HANGUL SYLLABLE *([^AEIOUWY]*)([AEIOUWY]+)([^AEIOUWY]*)$"
            if(length(p <- grep(re, x)))
                y[p] <- .u_char_from_name_Hangul(x[p])
            re <- "^CJK UNIFIED IDEOGRAPH ([0123456789ABCDEF]{4,6})$"
            p <- grep(re, x)
            if(length(p)) {
                u <- as.u_char(sub(re, "\\1", x[p]))
                ## Range check.
                p <- p[u %uin% UCD_Unicode_data_range$Range[1L : 4L]]
                y[p] <- u[p]
            }
        }
    } else {
        ## Currently, no support for Hangul syllable and CJK unified
        ## ideographs ...
        y <- with(UCD_Unicode_data_table,
                  Code[grepl(x, Name, ...)])
    }
        
    y
}

.u_char_from_name_Hangul <-
function(s)
{
    y <- rep.int(NA_integer_, length(s))
    re <- "^HANGUL SYLLABLE *([^AEIOUWY]*)([AEIOUWY]+)([^AEIOUWY]*)$"
    ok <- grepl(re, s)
    if(any(ok)) {
        s <- s[ok]
        L <- sub(re, "\\1", s)
        V <- sub(re, "\\2", s)
        T <- sub(re, "\\3", s)
        ## Seems we cannot subscript with an empty character string:
        y[ok] <- (Jamo_S_base
                  + (match(L, Jamo$L) - 1L) * Jamo_N_count
                  + (match(V, Jamo$V) - 1L) * Jamo_T_count
                  + (match(T, Jamo$T) - 1L)
                  )
    }
    as.u_char(y)
}

u_char_info <-
function(x)
{
    x <- as.u_char(x)
    ## Could make this more efficient eventually ...
    p <- match(x, UCD_Unicode_data_table$Code)
    y <- UCD_Unicode_data_table[p, ]
    row.names(y) <- NULL
    p <- which(is.na(p))
    if(length(p)) {
        r <- u_char_match(x[p], UCD_Unicode_data_range$Range, 0L)
        q <- p[r == UCD_Unicode_data_range_pos_Hangul]
        p <- p[r > 0L]
        y[p, ] <- UCD_Unicode_data_range[r, -1L]
        ## Not quite efficient ...
        y$Name[p] <- u_char_name(x[p])
        ## Call directly for efficiency ...
        y$Decomposition[q] <-
            format(.u_char_decomposition_Hangul(x[q]))
    }
    y
}

u_char_inspect <-
function(x)
{
    x <- as.u_char(x)
    data.frame(Code = x,
               Name = u_char_name(x),
               Char = intToUtf8(x, multiple = TRUE),
               stringsAsFactors = FALSE)
}

u_char_label <-
function(x)
{
    x <- as.u_char(x)
    ## Use the name by default.
    y <- u_char_name(x)
    ## Need to special case the ones with empty or missing names.
    p <- which(is.na(y) | (y == ""))
    if(length(p)) {
        ## Handle surrogate, private use and control characters.
        r <- u_char_match(x[p], 
                          c(UCD_Unicode_data_range$Range[-(1L : UCD_Unicode_data_range_pos_Hangul)],
                            "0000..001F", "007F..009F"),
                          0L)
        ## <FIXME>
        ## This hard-wires the positions of surrogate and private use
        ## character ranges.  Ideally, we would compute these (e.g.,
        ## from the rownames of UCD_Unicode_data_range).
        ## </FIXME>
        q <- p[(r > 0L) & (r <= 3L)]
        y[q] <- sprintf("surrogate-%04X", x[q])
        q <- p[(r > 3L) & (r <= 6L)]
        y[q] <- sprintf("private-use-%04X", x[q])
        q <- p[(r > 6L) & (r <= 8L)]
        y[q] <- sprintf("control-%04X", x[q])
        ## This leaves non-characters and reserved code points.
        q <- p[r == 0L]
        if(length(q)) {
            y[q] <-
                ifelse(u_char_property(x[q],
                                       "Noncharacter_Code_Point") ==
                       "Yes",
                       sprintf("noncharacter-%04X", x[q]),
                       sprintf("reserved-%04X", x[q]))
        }
    }
        
    y
}

u_char_name <-
function(x)
{
    x <- as.u_char(x)
    ## Could make this more efficient eventually ...
    p <- match(x, UCD_Unicode_data_table$Code)
    y <- UCD_Unicode_data_table$Name[p]
    p <- which(is.na(p))
    if(length(p)) {
        range <- UCD_Unicode_data_range$Range
        nomatch <- length(range) + 1L
        r <- u_char_match(x[p], range, nomatch)
        ph <- UCD_Unicode_data_range_pos_Hangul
        ind <- r < ph                   # CJK
        if(any(ind)) {
            pos <- p[ind]
            y[pos] <- .u_char_name_CJK(x[pos])
        }
        ind <- r == ph                  # Hangul
        if(any(ind)) {
            pos <- p[ind]
            y[pos] <- .u_char_name_Hangul(x[pos])
        }
        ind <- (r > ph) & (r < nomatch)
        if(any(ind)) {
            y[p[ind]] <- ""
        }
    }
    y

}

.u_char_name_CJK <-
function(x)
    sprintf("CJK UNIFIED IDEOGRAPH %04X", x)

.u_char_name_Hangul <-
function(x)
{
    s <- .u_char_decomposition_Hangul(x)
    L <- Jamo$L[sapply(s, `[`, 1L) - Jamo_L_base + 1L]
    V <- Jamo$V[sapply(s, `[`, 2L) - Jamo_V_base + 1L]
    T <- Jamo$T[sapply(s, `[`, 3L) - Jamo_T_base + 1L]
    T[is.na(T)] <- ""
    sprintf("HANGUL SYLLABLE %s%s%s", L, V, T)
}

u_char_property_names <-
function()
    unique(UCD_property_aliases_hash)

u_char_properties <-
function(x, which)
{
    if(nargs() == 0L) {
        ## As a convenience, list available properties if no arguments
        ## were given.
        return(sort(intersect(u_char_property_names(),
                              c(names(UCD_Unicode_data_table),
                                names(u_char_property_db),
                                names(u_char_property_getter_db)))))
    }
    
    x <- as.u_char(x)
    which <- .expand_property_aliases(as.character(which))
    all_property_names <- u_char_property_names()
    p <- match(which, all_property_names)
    if(any(ina <- is.na(p))) {
        warning(gettextf("Cannot map %s to Unicode properties",
                         paste(which[ina], collapse = ", ")),
                domain = NA)
        which <- which[!ina]
        p <- p[!ina]
    }

    ## Start with an empty data frame.
    y <- data.frame(x)[-1L]
    if(!any(duplicated(x)))
        row.names(y) <- x
    if(!length(which)) return(y)
    
    UCD_Unicode_data_table_names <- names(UCD_Unicode_data_table)
    p <- match(which, UCD_Unicode_data_table_names, 0L)
    if(any(p > 0L))
        y <- u_char_info(x)[p]
    ## Everything else would need to come from the property db
    p <- match(which, names(u_char_property_db), 0L)
    for(i in which(p > 0L))
        y[[which[i]]] <-
            .u_char_property_from_property_db(x, which[i])
    ## or the property getter db.
    p <- match(which, names(u_char_property_getter_db), 0L)
    for(i in which(p > 0L))
        y[[which[i]]] <- u_char_property_getter_db[[p[i]]](x)

    if(length(bad <- setdiff(which, names(y)))) {
        warning("Obtaining the following properties is not yet implemented:\n",
                strwrap(paste(bad, collapse = ", "),
                        indent = 4L, exdent = 4L))
        for(b in bad) y[[b]] <- rep.int(NA_character_, length(x))
    }

    y[match(which, names(y))]
}
    
u_char_property_getter_db <- list()

u_char_property_getter_db$Decomposition_Mapping <-
    u_char_decomposition

u_char_property <-
function(x, which)
{
    if(length(which <- as.character(which)) != 1L)
        stop(gettextf("Invalid '%s' argument.", "which"),
             domain = NA)
    u_char_properties(x, which)[[1L]]
}

.u_char_property_from_property_db <-
function(x, which)
{
    db <- u_char_property_db[[which]]
    y <- rep(db$Default, length(x))
    p <- u_char_match(x, db$Table[[1L]], 0L)
    y[p > 0L] <- db$Table[[2L]][p]
    .expand_property_value_aliases(y, which)
}

u_char_show_glyph <-
function(x)
{
    x <- as.u_char(x)
    ## For now we can only show one glyph at a time ...
    l <- length(x)
    if(l == 0L)
        message("Nothing to show.")
    else if(l > 1L) {
        message("Can only show one glyph.")
        x <- x[1L]
    }
    browseURL(sprintf("http://www.unicode.org/cgi-bin/refglyph?24-%04x",
                      x))
}
    
## See Section 3.12 of the Unicode Standard.

Jamo <-
    list(L =
         ## Initial (HANGUL CHOSEONG)
         c("G", "GG", "N", "D", "DD", "R", "M", "B", "BB", "S", "SS",
           "", "J", "JJ", "C", "K", "T", "P", "H"),
         V =
         ## Medial  (HANGUL JUNGSEONG)
         c("A", "AE", "YA", "YAE", "EO", "E", "YEO", "YE", "O", "WA",
           "WAE", "OE", "YO", "U", "WEO", "WE", "WI", "YU", "EU", "YI",
           "I"),
         T =
         ## Final   (HANGUL JONGSEONG)
         c("", "G", "GG", "GS", "N", "NJ", "NH", "D", "L", "LG", "LM",
           "LB", "LS", "LT", "LP", "LH", "M", "B", "BS", "S", "SS",
           "NG", "J", "C", "K", "T", "P", "H"))


Jamo_S_base <- hex("AC00")
Jamo_L_base <- hex("1100")
Jamo_V_base <- hex("1161")
Jamo_T_base <- hex("11A7")
Jamo_S_count <- 11172L
Jamo_L_count <- 19L
Jamo_V_count <- 21L
Jamo_T_count <- 28L
Jamo_N_count <- Jamo_V_count * Jamo_T_count

