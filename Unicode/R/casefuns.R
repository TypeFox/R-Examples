## Full Case Mappings.
## See Unicode 3.13.

## Use the map from SpecialCasing.txt, plus the maps from
## UnicodeData.txt, excluding any of the latter mappings that would
## conflict.  Any character that does not have a mapping in these files
## is considered to map to itself.

## <FIXME>
## For now, only use the unconditional maps in SpecialCasing.txt.
## Add arguments allowing for conditional maps eventually ...
## </FIXME>

u_to_lower_case <-
function(x)
    UseMethod("u_to_lower_case")

u_to_lower_case.u_char <-
function(x)
{
    y <- as.list(x)
    ## Maps from SpecialCasing.txt first.
    p <- match(x, UCD_special_casing_table$Code, 0L)
    p[p > 0L][nzchar(UCD_special_casing_table$Condition[p])] <- 0L
    y[p > 0L] <- UCD_special_casing_table$Lower[p]
    ## Maps from UnicodeData.txt not excluded by the above.
    q <- match(x, UCD_Unicode_data_table$Code, 0L)
    q[p > 0L] <- 0L
    r <- UCD_Unicode_data_table$Simple_Lowercase_Mapping[q]
    ind <- nzchar(r)
    y[q > 0L][ind] <- r[ind]

    as.u_char_seq(y)
}

u_to_lower_case.u_char_range <-
function(x)
    u_to_lower_case(as.u_char_seq(x))

u_to_lower_case.u_char_seq <-
function(x)
    as.u_char_seq(lapply(unclass(x),
                         function(e)
                         unlist(u_to_lower_case(e))))

u_to_lower_case.default <-
function(x)
{
    y <- lapply(x,
                function(s)
                intToUtf8(unlist(u_to_lower_case(.str_to_u_char(s)))))
    as.character(unlist(y))
}

u_to_upper_case <-
function(x)
    UseMethod("u_to_upper_case")

u_to_upper_case.u_char <-
function(x)
{
    y <- as.list(x)
    ## Maps from SpecialCasing.txt first.
    p <- match(x, UCD_special_casing_table$Code, 0L)
    p[p > 0L][nzchar(UCD_special_casing_table$Condition[p])] <- 0L
    y[p > 0L] <- UCD_special_casing_table$Upper[p]
    ## Maps from UnicodeData.txt not excluded by the above.
    q <- match(x, UCD_Unicode_data_table$Code, 0L)
    q[p > 0L] <- 0L
    r <- UCD_Unicode_data_table$Simple_Uppercase_Mapping[q]
    ind <- nzchar(r)
    y[q > 0L][ind] <- r[ind]

    as.u_char_seq(y)
}

u_to_upper_case.u_char_range <-
function(x)
    u_to_upper_case(as.u_char_seq(x))

u_to_upper_case.u_char_seq <-
function(x)
    as.u_char_seq(lapply(unclass(x),
                         function(e)
                         unlist(u_to_upper_case(e))))

u_to_upper_case.default <-
function(x)
{
    y <- lapply(x,
                function(s)
                intToUtf8(unlist(u_to_upper_case(.str_to_u_char(s)))))
    as.character(unlist(y))
}

u_to_title_case <-
function(x)
    UseMethod("u_to_title_case")

u_to_title_case.u_char <-
function(x)
{
    y <- as.list(x)
    ## Maps from SpecialCasing.txt first.
    p <- match(x, UCD_special_casing_table$Code, 0L)
    p[p > 0L][nzchar(UCD_special_casing_table$Condition[p])] <- 0L
    y[p > 0L] <- UCD_special_casing_table$Title[p]
    ## Maps from UnicodeData.txt not excluded by the above.
    q <- match(x, UCD_Unicode_data_table$Code, 0L)
    q[p > 0L] <- 0L
    r <- UCD_Unicode_data_table$Simple_Titlecase_Mapping[q]
    ind <- nzchar(r)
    y[q > 0L][ind] <- r[ind]

    as.u_char_seq(y)
}

## No other methods for now:
## For strings, we must find the word boundaries according to UAX #29
## "Unicode Text Segmentation", and then map characters following the
## word boundaries to their titlecase mapping, and the others to their
## lowercase mapping.


## Case folding.
## See Unicode 3.13.

## <FIXME>
## Add a status/mode argument eventually.
## For now, perform full case folding using mappings with status C and F.
u_case_fold <-
function(x)
    UseMethod("u_case_fold")
## </FIXME>

u_case_fold.u_char <-
function(x)
{
    y <- as.list(x)
    p <- match(x, UCD_case_folding_table$Code, 0L)
    p[(p > 0L) &
      is.na(match(UCD_case_folding_table$Status[p], c("C", "F")))] <- 0L
    y[p > 0L] <- UCD_case_folding_table$Mapping[p]

    as.u_char_seq(y)
}

u_case_fold.u_char_range <-
function(x)
    u_case_fold(as.u_char_seq(x))

u_case_fold.u_char_seq <-
function(x)
    as.u_char_seq(lapply(unclass(x),
                         function(e)
                         unlist(u_case_fold(e))))

u_case_fold.default <-
function(x)
{
    y <- lapply(x,
                function(s)
                intToUtf8(unlist(u_case_fold(.str_to_u_char(s)))))
    as.character(unlist(y))
}           
