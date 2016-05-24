
# internal function
evs.categories <- function(evs2008.lvl2, income="(all)", hhtype="(all)", relative=TRUE, ...) {
    dat <- evs2008.lvl2[(evs2008.lvl2$coicop2 != "15" & evs2008.lvl2$coicop2 != "00"),]
    income_lvls <- unique(dat$income)
    if(!income %in% income_lvls) stop("invalid 'income' argument, expected one of '", paste(income_lvls, collapse="', '"), "'.")
    hhtype_lvls <- unique(dat$hhtype)
    if(!hhtype %in% hhtype_lvls) stop("invalid 'hhtype' argument, expected one of '", paste(hhtype_lvls, collapse="', '"), "'.")
    dat <- dat[dat$income==income & dat$hhtype==hhtype,]
    if(relative) dat$value <- dat$value/sum(dat$value)
    res <- dat$value
    names(res) <- dat$coicop2de
    return(res)
}

# internal function
evs.elasticity <- function(evs2008.lvl2, categ="", xlab="", ylab="", main=NULL, ...) {
    #if(!require(beeswarm)) stop("could not load required package 'beeswarm'")
    dat <- evs2008.lvl2[(evs2008.lvl2$coicop2 != "15" & evs2008.lvl2$coicop2 != "00" & evs2008.lvl2$income != "(all)"),]
    cat_lvls <- unique(dat$coicop2)
    if (!categ %in% cat_lvls) stop("invalid 'categ' argument, expected one of '", paste(cat_lvls, collapse="', '"), "'.")
    
    income_lvls <- c("lt900", "900to1300", "1300to1500", "1500t2000", "2000t2600",
                    "2600t3600", "3600t5000", "5000t18000")
    dat$income <- factor(dat$income, levels=income_lvls)
    dat <- dat[dat$coicop2==categ,]
    if(is.null(main)) main <- dat[1, "coicop2de"]
    
    par(mar=c(2,2,3,0)+0.1)
    boxplot(value ~ income, data=dat, ylab=ylab, main=main, ylim=c(0, 1000), ...)
    
    return(dat)
}

#' Elasticity of Private Households
#'
#' This is a proof of concept for a use case of \code{datamart} to be used
#' for organising internal datasets as well as some calculations and graphics 
#' on these datasets.
#'
#' The example uses data from the Germany's Sample survey of income and expenditure (Einkommens- und Verbrauchsstichprobe, EVS).
#'
#' @seealso \code{\link{internalData}}, \code{\link{datamart}}
#'
#' @references 
#' Statistisches Bundesamt: Wirtschaftsrechnungen. Einkommens- und Verbrauchsstichprobe. Einnahmen und Ausgaben privater Haushalte. Fachserie 15 Heft 4.
#' @examples
#' xp <- expenditures()
#' queries(xp)
#' query(xp, "categories")
#' query(xp, "elasticity", categ="05")
#' @name expenditures
#' @export
expenditures <- function() datamart(
    internalData(name="evs2008.lvl2", package="datamart"),
    # resfunc(resource="elasticities", fun=evs.elasticities),
    resfunc(resource="categories", depends="evs2008.lvl2", fun=evs.categories),
    resfunc(resource="elasticity", depends="evs2008.lvl2", fun=evs.elasticity)
)
