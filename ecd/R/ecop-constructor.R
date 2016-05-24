#' Constructor of ecop class by read conf for option sample data
#'
#' Read conf for option sample data and fitting parameters
#'
#' @param key    character. The top-level key in conf
#' @param conf_file file name fof symbol config, default to conf/ecld-fit-conf.yml
#' @param conf_data optionally feed config through a list.
#'                  If this is not null, this takes priority and \code{conf_file} will be ignored.
#' @param ecop an ecop object with conf
#' @param df dataframe of a single closing date and time to maturity
#' @param otype option type
#'
#' @return the ecop object
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecop.from_symbol_conf
#' @export ecop.read_symbol_conf
#'
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#'     conf <- ecop.read_symbol_conf("spx2_1d")
#'     op <- ecop.from_symbol_conf("spx2_1d")
#' }
### <======================================================================>
"ecop.from_symbol_conf" <- function(key,
                                    conf_file = "conf/ecop-fit-conf.yml",
                                    conf_data=NULL)
{
    # read the conf data
    conf <- NULL
    if (! is.null(conf_data) & class(conf_data)=="list") {
        conf <- conf_data
    } else {
        conf <- ecop.read_symbol_conf(key, conf_file)
    }
    if (is.null(conf$symbol)) {
        stop("Failed to locate symbol in conf")
    }
    # basic info
    if (!("int_rate" %in% names(conf))) {
        conf$int_rate <- 0.0
    }
    if (!("div_yield" %in% names(conf))) {
        conf$div_yield <- 0.0
    }
    if (!("ttm" %in% names(conf))) {
        conf$ttm <- conf$days/365
    }
    # fill up empty list for missing call or put
    if (!("call" %in% names(conf))) {
        conf$call <- list()
    }
    if (!("put" %in% names(conf))) {
        conf$put <- list()
    }
    #
    call_conf <- conf$call
    put_conf <- conf$put # if anything missing, always get it from call conf
    
    if (!("ecld" %in% names(put_conf)) & ("ecld" %in% names(call_conf))) {
        put_conf$ecld <- call_conf$ecld
    }
    
    # first pass of ecop object
    call <- match.call()
    op = new("ecop", call = call,
        conf = conf,
        key = key,
        symbol = conf$symbol,
        datadate = as.Date(conf$datadate),
        days = conf$days,
        ttm = conf$ttm2[1] / conf$ttm2[2],
        int_rate = conf$int_rate,
        div_yield = conf$div_yield,
        put_conf = put_conf,
        call_conf = call_conf,
        put_data = new("ecop.opt"),
        call_data = new("ecop.opt")
    )
    # option data
    Date = days = NULL # Simply to avoid R CMD check complaining
    df <- ecop.read_csv_by_symbol(op@symbol)
    df2 <- subset(df, Date == op@datadate & days == op@days )
    df2$diff <- df2$STRK_PRC - df2$UNDL_PRC
    
    if (length(names(call_conf)) > 0) {
        op@call_data <- ecop.build_opt(op, df2, otype="c")
    }
    if (length(names(call_conf)) > 0) {
        op@put_data <- ecop.build_opt(op, df2, otype="p")
    }
    # done
    invisible(op)
}
### <---------------------------------------------------------------------->
#' @rdname ecop.from_symbol_conf
"ecop.read_symbol_conf" <- function(key, conf_file = "conf/ecop-fit-conf.yml")
{
    # read the conf file
    if (! file.exists(conf_file)) {
        stop(paste("conf_file does not exist:", conf_file))
    }
    conf_all <- yaml.load_file(conf_file)
    conf1 <- conf_all[[key]]
    conf1
}
### <---------------------------------------------------------------------->
#' @rdname ecop.from_symbol_conf
# this is an internal utility
"ecop.build_opt" <- function(ecop, df, otype)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    conf <- NULL
    if (otype=="c") conf <- ecop@call_conf
    if (otype=="p") conf <- ecop@put_conf
    
    PC = NULL # Simply to avoid R CMD check complaining
    df2 <- subset(df, diff >= conf$range.from & diff <= conf$range.to)
    df2 <- subset(df2, tolower(PC)==otype)
    K <- df2$STRK_PRC
    S <- S_raw <- head(unique(df2$UNDL_PRC),1)
    if ("S_override" %in% names(conf)) {
        S <- conf$S_override
    }
    
    r <- ecop@int_rate
    T <- ecop@ttm
    
    ld <- ecld(sigma=1) # meaningless default
    if ("ecld" %in% names(conf)) {
        ldc <- conf$ecld
        is.sged <- if ("is_sged" %in% names(ldc)) ifelse(ldc$is_sged > 0, TRUE, FALSE) else FALSE
        ld0 <- ecld(lambda=ldc$lambda, sigma=ldc$sigma, beta=ldc$beta, is.sged=is.sged)
        mu_D <- ecld.mu_D(ld0)
        mu <- mu_D + ldc$mu_plus
        ld <- ecld(lambda=ldc$lambda, sigma=ldc$sigma, beta=ldc$beta, mu=mu, is.sged=is.sged)
    }
    if ("ecdr" %in% names(conf)) {
        ldc <- conf$ecdr
        ld0 <- ecd.polar(R=ldc$R, theta=ldc$degree/180*pi, sigma=ldc$sigma*ecd.mp1, beta=ldc$beta)
        mu_D <- ecd.mu_D(ld0)
        mu <- mu_D + ldc$mu_plus
        ld <- ecd.polar(R=ldc$R, theta=ldc$degree/180*pi, sigma=ldc$sigma*ecd.mp1, beta=ldc$beta, mu=mu)
    }
    
    if (!("momentum" %in% names(conf))) {
        conf$momentum <- 0.0
    }
    if (!("epsilon" %in% names(conf))) {
        conf$epsilon <- 0.0
    }
    if (!("k_cusp" %in% names(conf))) {
        conf$k_cusp <- 0.0
    }
    
    call <- match.call()
    opt <- new("ecop.opt", call = call,
        otype = otype,
        range.from = conf$range.from,
        range.to = conf$range.to,
        momentum = conf$momentum,
        epsilon = conf$epsilon,
        k_cusp = conf$k_cusp,
        ecldOrEcd = ld,
        S = S,
        S_raw = S_raw,
        strike = K,
        k = log(K/S) - r*T,
        V_last = df2$LAST,
        V_bid = df2$L_BID,
        V_ask = df2$L_ASK,
        V = df2$L_ASK/2 + df2$L_BID/2
    )
    opt
}
### <---------------------------------------------------------------------->
