rsmformula <- function(design, response=NULL, factor.names=NULL, 
        use.blockvar = TRUE, degree=2, coded=TRUE, ...){
      ##dnam <- deparse(substitute(design))

      ## avoid NOTE about no global function definition
      ##if (!is.loaded("FrF2")) iscube <- function(design) NULL
      ## does not work


      if (!"design" %in% class(design)) stop("rsmformula is applicable to class design objects only.")
      di <- design.info(design)
      if (coded & is.null(di$coding))
          stop("rsmformula with coded values is applicable to designs with coding information only.")
      if (is.null(di$response.names)) stop("rsmformula needs at least one response in the design")
      if (is.null(di$quantitative) | !all(di$quantitative)) stop("For rsmformula, all design factors must be quantitative.")
      if(!degree %in% c(1,1.5,2)) stop("rsmformula supports degrees 1, 1.5 and 2 only")
      if (!is.null(response))
         if (length(response)>1) stop("rsmformula can only handle one response at a time")
      if (is.null(factor.names)) factor.names <- names(di$factor.names)
      if (length(factor.names)<2) stop("rsmformula needs at least two experimental factors")
      if (!is.logical(use.blockvar)) stop("use.blockvar must be a logical")

      ## check available responses
      respnam <-di$response.names
      respnamOK <- intersect(colnames(design),respnam)
      if (is.null(respnamOK) | length(respnamOK)==0)
          stop("For formula.design, the design requires at least one response to be available.")
      respnamOK <- respnamOK[which(sapply(design[,respnamOK], function(obj) all(!is.na(obj))))]
      if (length(respnamOK)==0) stop("the design does not contain any response variable with complete observations")
      respposOK <- which(di$response.names %in% respnamOK)

      ## check response given by user
      if (!is.null(response)){
         if (!(is.character(response) | is.numeric(response)))
              stop("response must be a character string of the response name or a position number")
         if (is.numeric(response)) {
               if (response < 1 | response > length(respnam) | !response==round(response))
                     stop("if numeric, response must be an integer from 1 to ", length(respnam))
               response <- respnam[response]
               ## now response is character
            }
         if (is.character(response)){
               if (!response %in% colnames(design))
                     stop("response is not a column of design")
               if (!response %in% respnam)
                     stop("response has not been declared a response variable")
               if (!response %in% respnamOK)
                     stop("response has missing values, which precludes default analysis of the design")
            }
         }
      else response <- respnamOK[1]
      ## else: no response given by user
      
      ## check factor.names 
      if (!is.character(factor.names))
              stop("factor.names must be a character vector")
      if (!all(factor.names %in% names(di$factor.names)))
              stop("invalid name(s) included in factor.names")
      which.facs <- which(names(di$factor.names) %in% factor.names)

      ## create RHS formula
      vars <- paste(names(di$factor.names[which.facs]),collapse=",")
      if (coded & !is.null(di$coding)) vars <- paste(names(di$coding)[which.facs],collapse=",")
      if (degree==1) rhstxt <- paste("FO(", vars,")", sep="") 
      else if (degree==1.5) rhstxt <- paste("FO(", vars, ") + TWI(", vars,")", sep="")
           else if (degree==2) rhstxt <- paste("FO(", vars, ") + TWI(", vars,") + PQ(", vars,")", sep="")
      ## add block if requested
      if (use.blockvar) if (!is.null(di$block.name)) rhstxt <- paste(di$block.name, rhstxt, sep="+")
      formula(paste(response, rhstxt, sep="~"))
}