grofit.control <-
function(
neg.nan.act       = FALSE,
clean.bootstrap   = TRUE ,
suppress.messages = FALSE,
fit.opt     = "b",
log.x.gc    = FALSE,
log.y.gc    = FALSE,
interactive = TRUE,
nboot.gc    = 0,
smooth.gc   = NULL,
model.type  = c("logistic", "richards", "gompertz", "gompertz.exp"),
have.atleast   = 6,
parameter      = 9,
smooth.dr      = NULL,
log.x.dr       = FALSE,
log.y.dr       = FALSE,
nboot.dr      = 0)
{
if ((is.character(fit.opt)==FALSE)|(length(fit.opt)!=1))
    stop("value of fit.opt must be character and of one element")

if (is.character(model.type)==FALSE)
    stop("value of model.type must be character")

if ((is.logical(neg.nan.act)==FALSE)|(length(neg.nan.act)!=1))
    stop("value of neg.nan.act must be logical and of one element")

if ((is.logical(clean.bootstrap)==FALSE)|(length(clean.bootstrap)!=1))
    stop("value of clean.bootstrap must be logical and of one element")

if ((is.logical(suppress.messages)==FALSE)|(length(suppress.messages)!=1))
    stop("value of suppress.messages must be logical and of one element")

if ((is.logical(log.x.gc)==FALSE)|(length(log.x.gc)!=1))
    stop("value of log.x.gc must be logical and of one element")

if ((is.logical(log.y.gc)==FALSE)|(length(log.y.gc)!=1))
    stop("value of log.y.gc must be logical and of one element")

if ((is.logical(interactive)==FALSE)|(length(interactive)!=1))
    stop("value of interactive must be logical and of one element")

if ((is.logical(log.x.dr)==FALSE)|(length(log.x.dr)!=1))
    stop("value of log.x.dr must be logical and of one element")

if ((is.logical(log.y.dr)==FALSE)|(length(log.y.dr)!=1))
    stop("value of log.y.dr must be logical and of one element")

if ((is.numeric(nboot.gc)==FALSE)|(length(nboot.gc)!=1)|(nboot.gc<0))
    stop("value of nboot.gc must be numeric (>=0) and of one element")

if ((is.numeric(have.atleast)==FALSE)|(length(have.atleast)!=1)|(have.atleast<6))
    stop("value of have.atleast must be numeric (>=6) and of one element")

if ((is.numeric(parameter)==FALSE)|(length(parameter)!=1))
    stop("value of parameter must be numeric and of one element")

if ((is.numeric(nboot.dr)==FALSE)|(length(nboot.dr)!=1)|(nboot.dr<0))
    stop("value of nboot.dr must be numeric (>=0) and of one element")

if (((is.numeric(smooth.gc)==FALSE) && (is.null(smooth.gc)==FALSE)))
    stop("value of smooth.gc must be numeric or NULL")

if (((is.numeric(smooth.dr)==FALSE) && (is.null(smooth.dr)==FALSE)))
    stop("value of smooth.dr must be numeric or NULL")


grofit.control <- list(neg.nan.act=neg.nan.act, clean.bootstrap=clean.bootstrap, suppress.messages=suppress.messages,
                       fit.opt=fit.opt, log.x.gc=log.x.gc, log.y.gc=log.y.gc, interactive=interactive,
                       nboot.gc=round(nboot.gc), smooth.gc=smooth.gc, smooth.dr=smooth.dr,
                       have.atleast=round(have.atleast), parameter=round(parameter), log.x.dr=log.x.dr, log.y.dr=log.y.dr,
                       nboot.dr=round(nboot.dr), model.type=model.type)
class(grofit.control) <- "grofit.control"
grofit.control
}

