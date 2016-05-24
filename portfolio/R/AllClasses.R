################################################################################
##
## $Id: AllClasses.R 1311 2008-10-31 17:38:03Z enos $
##
## Class definitions for the portfolio package.
##
################################################################################

## Load hook for methods

setOldClass(c("Date"))

## The size slot of the portfolioBasic class requires either a
## character or numeric vector.

setClassUnion("characterOrNumeric", c("character","numeric"))
setClassUnion("optionalList", c("list", "NULL"))
setClassUnion("optionalFormula", c("formula", "NULL"))

## The portfolioBasic class contains the basic set of information
## required for working with portfolios with only weight data.

setClass("portfolioBasic",
         representation(
                        name          = "character",

                        ## The instant to which this portfolio
                        ## pertains.  Currently this slot can contain
                        ## anything -- is there any reason to restrict
                        ## it?
                        
                        instant       = "ANY",

                        ## The data slot can contain anything, but
                        ## must contain a column named 'id'.  Some
                        ## portfolio methods will require that other
                        ## columns be present in this data frame, and
                        ## will fail if they are missing.
                        
                        data          = "data.frame",

                        ## The following variable specs are character
                        ## vectors of length 1.  There can be only one
                        ## 'ret.var', for instance.

                        id.var        = "character",
                        symbol.var    = "character",
                        in.var        = "character",
                        weight.var    = "character",
                        ret.var       = "character",

                        type          = "character",
                        size          = "characterOrNumeric",
                        weight.range  = "numeric",
                        sides         = "character",
                        
                        ## Weights storage.  Contains data frame with
                        ## names c("id","weight").
                        
                        weights       = "data.frame"
                        ),

         ## in.var and weight.var are zero-length vectors to avoid a
         ## class union and to ease delegation to the 'weight'
         ## function.  It would be incorrect to default to some
         ## potentially missing (or unwanted) variable.
         
         prototype = prototype(
           name          = "Unnamed portfolio",
           instant       = NULL,
           data          = data.frame(id = I(character(0))),

           id.var        = "id",
           symbol.var    = character(0),
           in.var        = character(0),
           weight.var    = character(0),
           ret.var       = character(0),
           
           type          = "equal",
           size          = "quintile",
           weight.range  = c(0, Inf),
           sides         = c("long","short"),
           
           weights       = data.frame(id = I(character(0)), weight = numeric(0))
           ),

         validity = function(object){

           if(length(object@id.var) != 1){
             return("id.var slot must be of length 1.")
           }
           
           if(!all(c("id","weight") %in% names(object@weights))){
             return(paste("Data frame in slot 'weights' must contain",
                          "columns 'id' and 'weight'"))
           }
           if(!is.character(object@weights$id)){
             return("Invalid class for id column in weights slot: must be character")
           }
           if(any(is.na(object@weights$id))){
             return("NA's in column 'id' in weights slot not allowed")
           }
           if(!is.numeric(object@weights$weight)){
             return("Invalid class for weight column in weights slot: must be numeric")
           }
           if(any(duplicated(object@data$id))){
             return("Cannot have duplicate id's in the data slot")
           }
           if(any(duplicated(object@weights$id))){
             return("Cannot have duplicate id's in the weights slot")
           }

           ## Check for illegal column names in the data slot

           illegal.cols <- c("weight")
           if(any(illegal.cols %in% names(object@data))){
             illegal.present <- illegal.cols[illegal.cols %in% names(object@data)]
             return(paste("The following columns are not allowed in the data slot:",
                          paste(illegal.present, collapse = ",")))
           }
           
           w <- object@weights
           w <- w[order(w$id),]

           ## I'm torn on the issue of whether all weights/shares
           ## entries need to have corresponding data slot entries.
           ## I'm back to looser constraints because of my tendency to
           ## update the entire data slot without being terribly
           ## careful about the securities in it.

           ## As a compromise, warn for now.

           if(!all(w$id %in% object@data$id)){
             id.viol <- w$id[which(!w$id %in% object@data$id)]
             
             warning(paste("The securities with the following id's have",
                           "weights but no data in the data slot:", id.viol))
           }
           
           return(TRUE)
         }
         )

## The portfolio class extends the portfolioBasic class for working
## with price and share information.

setClass("portfolio",
         representation(
                        equity       = "numeric",
                        weight.style = "character",
                        file         = "character",
                        price.var    = "character",
                        shares       = "data.frame"
                        ),
         prototype = prototype(
           equity       = numeric(0),
           weight.style = "sides.separate",
           file         = "none",
           price.var    = "price.usd",
           
           shares = data.frame(id = I(character(0)), shares = numeric(0))
           ),
         contains = "portfolioBasic",

         validity = function(object){
           if(!all(c("id","shares") %in% names(object@shares))){
             return(paste("Data frame in slot 'shares' must contain",
                          "columns 'id' and 'shares'"))
           }
           if(!is.character(object@shares$id)){
             return("Invalid class for id column in shares slot: must be character")
           }
           if(any(is.na(object@shares$id))){
             return("NA's in column 'id' in shares slot not allowed")
           }
           if(!is.numeric(object@shares$shares)){
             return("Invalid class for shares column in shares slot: must be numeric")
           }
           if(any(duplicated(object@shares$id))){
             return("Cannot have duplicate id's in the weights slot")
           }

           ## Some controls on weight.style and equity:

           if(length(object@weight.style) != 1){
             return("Slot weight.style must have length 1")
           }
           if(!object@weight.style %in% c("sides.separate", "long.tmv",
                                          "short.tmv", "reference.equity")){
             return("Invalid value in weight.style slot.  Please see class?portfolio.")
           }
           if(object@weight.style %in% "reference.equity" &&
              length(object@equity) != 1){
             return("equity slot must have length 1 for weight.type 'reference.equity'")
           }

           ## NA shares are not allowed.

           if(any(is.na(object@shares$shares))){
             return("NA shares are not allowed")
           }

           ## 0 shares not allowed (we may want to loosen this
           ## constraint later).

           if(any(object@shares$shares == 0)){
             return("0 shares not allowed")
           }
           
           ## Valid portfolio objects must have the same set of
           ## identifiers in the shares and weights data.frame's.

           w <- object@weights
           w <- w[order(w$id),]

           s <- object@shares
           s <- s[order(s$id),]

           ## Furthermore, the set of securities in the weights and
           ## shares must be all.equal, up to their character
           ## representations.
           
           return(all.equal(as.character(s$id), as.character(w$id)))
         }
         )

setClassUnion("portfolioBasicOrNull", c("portfolioBasic", "NULL"))
setClassUnion("portfolioOrNull", c("portfolio", "NULL"))

setClass("exposure",
         representation(
                        data  = "list"
                        ),
         prototype = prototype(
           data = list()
           )
         )

## Here we have ret, which is a unitless measure of return, and
## profit, which is a monetary amount reflecting profit and loss.  In
## retrospect, 'ret.detail' should have been called 'detail'.

setClass("performance",
         representation(
                        ret            = "numeric",
                        profit         = "numeric",
                        missing.price  = "numeric",
                        missing.return = "numeric",
                        ret.detail     = "data.frame",
                        t.plus.one     = "portfolioBasicOrNull"
                        ),
         prototype = prototype(
           ret            = 0,
           profit         = 0,
           missing.price  = 0,
           missing.return = 0,
           ret.detail     = data.frame(),
           t.plus.one     = NULL
           )
         )

setClass("contribution",
         representation(
                        data  = "list"
                        ),
         prototype = prototype(
           data = list()
           )
         )

## An object of the class "trades" contains a data frame with columns
## "id", "side", and "shares" describing a simple list of trades to be
## performed.

setClass("trades",

         representation(trades = "data.frame"),

         prototype(trades = data.frame(
                     id = character(),
                     side = character(),
                     shares = numeric())),

         validity = function(object){

           if(!isTRUE(all(c("id","side","shares") %in% names(object@trades)))){
             return("Columns \"id\", \"side\", and \"shares\" required.")
           }

           if(!isTRUE(all(!is.na(object@trades$id))) ||
              !isTRUE(all(!is.na(object@trades$side))) ||
              !isTRUE(all(!is.na(object@trades$shares)))){
             return("No NAs allowed in id, side, and shares columns.")
           }

           if(!isTRUE(all(object@trades$side %in% c("B","S","C","X")))){
             return("Sides must be one of \"B\", \"S\", \"C\", or \"X\"")
           }

           if(!isTRUE(all(object@trades$shares > 0))){
             return("Shares must be greater than 0.")
           }
           
           if(!is.numeric(object@trades$shares)){
             return("Values in the shares column must be numeric.")
           }

           if(!isTRUE(all(!duplicated(paste(object@trades$id,
                                            ifelse(object@trades$side %in% c("B","S"),
                                                   "long","short")))))){
             return("Only one trade per id per side allowed.")
           }

           return(TRUE)
         }
         
         )

setClass("tradelist",
         representation(
                        type          = "character",
                        
                        ## Intermediate steps of the process
                        ## are stored in the following slots:

                        candidates    = "data.frame",
                        ranks         = "data.frame",
                        chunks        = "data.frame",
                        swaps         = "data.frame",
                        swaps.actual  = "data.frame",
                        chunks.actual = "data.frame",

                        ## The actual slot is the data frame
                        ## from with an xml tradelist is created.

                        actual        = "data.frame",

                        ## The final slot is a simpler version of
                        ## actual, and only contains id, side, and
                        ## shares.

                        final         = "trades",
                        
                        ## Allows the user to specify the names of the
                        ## columns in "data" containing necessary
                        ## information for tradelist construction
                        
                        id.var        = "character",
                        price.var     = "character",
                        
                        ## The "sorts" slot contains named list of
                        ## sorts. The value paired with each sort is
                        ## the weight to be applied to the ranks of
                        ## the stocks in that sort. A weight of 1/10
                        ## means that 1 trade will appear for every 10
                        ## trades in a sort with weight 1.

                        sorts         = "optionalList",

                        ## The set of trades that pass
                        ## the sort criteria is stored in the
                        ## 'rank.sorts' slot after those trades are
                        ## computed.

                        rank.sorts    = "list",
                        
                        regions       = "character",
                        chunk.usd     = "numeric",
                        trade.usd.min = "numeric",

                        ## The 'unrestricted' flag directs processing
                        ## of this tradelist to ignore all possible
                        ## restrictions.  Trades considered are all
                        ## trades to proceed from the orig to target
                        ## portfolios.

                        unrestricted  = "logical",
                        
                        ## The 'restrictions' data frame specifies
                        ## which types of trades are not permitted,
                        ## and may be set by the client.

                        restrictions  = "data.frame",

                        ## The 'restricted' data frame contains
                        ## candidate trades removed due to a row in
                        ## the restrictions data frame, including
                        ## reason for removal.
                        
                        restricted     = "data.frame",

                        to.equity      = "logical",
                        turnover       = "numeric",
                        tca            = "character",
                        rank.gain.min  = "numeric",

                        
                        ## Here we store off some important
                        ## information about the target and original
                        ## portfolio, necessary for creating the
                        ## tradelist.

                        ## Why not allow a different target equity for
                        ## the long and short sides?  I believe the
                        ## correct answer is to provide for a long and
                        ## short reference equity in the portfolio
                        ## class.
                        
                        target.equity  = "numeric",
                        mv.long.orig   = "numeric",
                        mv.short.orig  = "numeric",

                        ## All other data goes here.
                        
                        data           = "data.frame",

                        ## Should methods on this object be verbose?

                        verbose        = "logical"
                        
                        ),

         prototype = prototype(
           type           = "ranks",
           candidates     = data.frame(),
           ranks          = data.frame(),
           chunks         = data.frame(),
           swaps          = data.frame(),
           swaps.actual   = data.frame(),
           chunks.actual  = data.frame(),
           actual         = data.frame(),
           final          = new("trades"),
           
           id.var         = "id",
           price.var      = "price.usd",
           
           sorts          = list(default.sort = 1),
           rank.sorts     = list(),
           
           regions        = character(0),

           unrestricted   = FALSE,
           restrictions   = data.frame(),
           restricted     = data.frame(),
           chunk.usd      = 10000,
           trade.usd.min  = 0,
           to.equity      = TRUE,
           turnover       = 0,
           tca            = c("volume"),
           rank.gain.min  = -Inf,
           data           = data.frame(),
           verbose        = FALSE
           ),
         validity = function(object){

           if(!object@type %in% c("ranks","all")){
             return(paste("Invalid type", object@type))
           }
           
           reserved <- c("orig", "target", "side", "shares", "mv")

           if(isTRUE(all.equal(object@type, c("ranks")))){

             reserved <- c(reserved, c("rank", "rank.t", "tca.rank",
                                       "chunk.shares", "chunk.mv", "chunk", "id.enter",
                                       "orig.enter", "target.enter", "side.enter",
                                       "shares.enter", "mv.enter", "rank.t.enter",
                                       "tca.rank.enter", "chunk.shares.enter", "chunk.mv.enter",
                                       "chunk.enter", "id.exit", "orig.exit", "target.exit",
                                       "side.exit", "shares.exit", "mv.exit", "rank.t.exit",
                                       "tca.rank.exit", "chunk.shares.exit", "chunk.mv.exit",
                                       "chunk.exit", "rank.gain"))
           }

           if(any(reserved %in% names(object@data))){
             return("The following column names are not allowed in \"data\":",
                    reserved[reserved %in% names(object@data)], "\n")
           }

           if(!object@id.var %in% names(object@data)){
             return("Column specified by \"id.var\" does not exist in \"data\".\n")
           }

           if(!c(object@price.var %in% names(object@data))){
             return("Column specified by \"price.var\" does not exist in \"data\".\n")
           }

           ## Currently, volume is only required for a tradelist of
           ## type ranks.

           if(isTRUE(all.equal(object@type, c("ranks")))){
             if(!c("volume" %in% names(object@data))){
               return("\"data\" requires a \"volume\" column.\n")
             }
           }
           TRUE
         }

         )

## Should matchedPortfolio extend portfolioBasic?

setClassUnion("formulaOrNull", c("formula", "NULL"))
setClass("matchedPortfolio",

         representation(formula  = "formulaOrNull",
                        method   = "character",
                        original = "portfolioBasic",
                        omitted.treatment  = "numeric",
                        omitted.control  = "numeric",
                        matches  = "matrix"
                        ),
                        
         prototype = prototype(
           formula  = NULL,
           method   = "random",
           original = new("portfolioBasic"),
           omitted.treatment = 0,
           omitted.control = 0,
           matches  = matrix()
           ),
         
         validity = function(object){
           
           ## The object is empty if the formula is null

           if(!is.null(object@formula)){

             covariates <- all.vars(getCovariateFormula(object@formula))
           
             ## All terms of the formula must have columns in the
             ## 'data' slot of 'original'

             if(!isTRUE(all(covariates %in%
                            names(object@original@data)))){

               return(paste("'data' slot in 'original' does not contain columns",
                            "for all variables in the formula."))

             }
             
             ## Make sure ret.var points to something.

             if(!isTRUE(object@original@ret.var %in%
                        names(object@original@data))){
               
               return(paste("'data' slot in 'original' does not contain a column for",
                            "the 'ret.var'"))

             }

           }

           ## Matches and the original portfolio should have the same number of stocks.

           if(!isTRUE(all.equal(nrow(object@original@weights),
                                dim(object@matches)[1]))){
             
             return(paste("Number of stocks in portfolio differs between 'original'",
                          "and 'matches'"))

           }

           return(TRUE)

         }
         )

setClass("matchedPortfolioCollection",
         representation(data = "list"),
         prototype = prototype(
           data = list()
           ),

         validity = function(object){
           if(!all(sapply(object@data, class) == "matchedPortfolio")){
             return("All elements must be of class 'matchedPortfolio'")
           }
           return(TRUE)
         }
         )
