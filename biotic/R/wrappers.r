#' Calculate BMWP indices for invertebrate samples
#' @description Calculates BMWP, ASPT and N-taxa index values for invertebrate samples.
#' @param df A dataframe containing list of taxa in first column, followed by
#'  columns of abundances with sample names in header row.
#' @param type Indicates format of data. Options are "num" for numeric data,
#' "log" for integer log abundance categories (1-5) or "alpha" for alphabetic
#' abundance categories (A-E). Default value is "num".
#' @return A data frame consisting of columns of index values with samples in
#' rows.
#' @export calcBMWP
#' @examples
#'
#' # calculate the BMWP indices for the River Almond dataset
#' # 'type' not specified as data are numeric abundances
#'
#' calcBMWP(almond)

calcBMWP<-function(df, type="num"){
  calcindex(df, index="BMWP", type=type)
}

#' Calculate Whalley revised BMWP indices for invertebrate samples
#' @description Calculates Whalley revised BMWP, ASPT and N-taxa indices
#' for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of columns of index values with samples
#' in rows.
#' @export calcWhalley
#' @examples
#'
#' # calculate the Whalley revised BMWP indices for the Green Burn dataset
#' # data are numeric log abundance categories, so type is "log"
#'
#' calcWhalley(greenburn, "log")

calcWhalley<-function(df, type="num"){
  calcindex(df, index="Whalley", type=type)
}

#' Calculate Whalley 'Riffle' habitat-specific BMWP indices for invertebrate
#' samples
#' @description Calculates Whalley riffle-specific BMWP, ASPT and N-taxa
#' indices for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of columns of index values with samples
#' in rows.
#' @export calcRiffle
#' @examples
#'
#' # calculate the Whalley Riffle BMWP indices for the Braid Burn dataset
#' # data are alphabetic log abundance categories, so type is "alpha"
#'
#' calcRiffle(braidburn, "alpha")


calcRiffle<-function(df, type="num"){
  calcindex(df, index="Riffle", type=type)
}

#' Calculate Whalley 'Pool' habitat-specific BMWP indices for invertebrate
#' samples
#' @description Calculates Whalley pool-specific BMWP, ASPT and N-taxa indices
#' for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of columns of index values with samples
#' in rows.
#' @export calcPool
#' @examples
#'
#' # calculate the Whalley Pool BMWP indices for the Green Burn dataset
#' # data are numeric log abundance categories, so type is "log"
#'
#' calcPool(greenburn, "log")

calcPool<-function(df, type="num"){
  calcindex(df, index="Pool", type=type)
}

#' Calculate Whalley 'Riffle/Pool' habitat-specific BMWP indices
#' @description Calculates Whalley riffle/pool-specific BMWP, ASPT and N-taxa
#' indices for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of columns of index values with samples
#' in rows.
#' @export calcRiffPool
#' @examples
#'
#' # calculate the Whalley Riffle/Pool BMWP indices for the Braid Burn dataset
#' # data are alphabetic log abundance categories, so type is "alpha"
#'
#' calcRiffPool(braidburn, "alpha")

calcRiffPool<-function(df, type="num"){
  calcindex(df, index="RiffPool", type=type)
}

#' Calculate presence-only WHPT indices
#' @description Calculates WHPT presence-only ASPT and N-taxa indices
#' for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of columns of index values with samples
#' in rows.
#' @export calcWHPT
#' @examples
#'
#' # calculate the WHPT presence-only indices for the Braid Burn dataset
#' # data are alphabetic log abundance categories, so type is "alpha"
#'
#' calcWHPT(braidburn, "alpha")

calcWHPT<-function(df, type="num"){
  calcindex(df, index="WHPT", type=type)
}

#' Calculate abundance-weighted WHPT indices
#' @description Calculates WHPT abundance-weighted ASPT and N-taxa indices
#' for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of columns of index values with samples
#' in rows.
#' @export calcWHPT_AB
#' @examples
#'
#' # calculate the WHPT abundance-weighted indices for the River Almond dataset
#' # data are numeric abundances, so type is "num" (can be omitted)
#'
#' calcWHPT_AB(almond, "num")

calcWHPT_AB<-function(df, type="num"){
  calcindex(df, index="WHPT_AB", type=type)
}

#' Calculate LIFE index
#' @description Calculates LIFE index values for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of a column of index values with samples
#' in rows.
#' @export calcLIFE
#' @examples
#'
#' # calculate the LIFE index for the River Almond dataset
#' # data are numeric abundances, so type can be omitted ("num" is default)
#'
#' calcLIFE(almond)

calcLIFE<-function(df, type="num"){
  calcindex(df, index="LIFE", type=type)
}


#' Calculate PSI index
#' @description Calculates PSI index for invertebrate samples.
#' for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of a column of index values with samples
#' in rows.
#' @export calcPSI
#' @examples
#'
#' # calculate the PSI index for the Braid Burn dataset
#' # data are alphabetic log abundance categories, so type is "alpha"
#'
#' calcPSI(braidburn, "alpha")

calcPSI<-function(df, type="num"){
  calcindex(df, index="PSI", type=type)
}

#' Calculate AWIC index
#' @description Calculates AWIC index for invertebrate samples.
#' @inheritParams calcBMWP
#' @return A data frame consisting of a column of index values with samples
#' in rows.
#' @export calcAWIC
#' @examples
#'
#' # calculate the AWIC index for the Green Burn dataset
#' # data are numeric log abundance categories, so type is "log"
#'
#' calcAWIC(greenburn, "log")

calcAWIC<-function(df, type="num"){
  calcindex(df, index="AWIC", type=type)
}
