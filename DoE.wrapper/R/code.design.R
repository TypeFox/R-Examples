.parse.coding <- function (form) 
{
    if (!inherits(form, "formula")) 
        stop("Coding formulas must be of class \"formula\"")
    if (length(form) < 3) 
        stop("Formula lacks a left-hand-side")
    nm = all.vars(form)
    if (length(nm) < 2) 
        stop(paste("Error in coding formula:", .form2str(form), 
            "\nCoded and uncoded names must differ"))
    names(nm) = c("coded", "orig")
    rhs = as.character(form)[3]
    a = eval(parse(text = sub(nm[2], "0", rhs)))
    b = eval(parse(text = sub(nm[2], "1", rhs)))
    d = 1/(b - a)
    list(names = nm, const = c(center = signif(-a * d, 4), divisor = signif(d, 
        4)))
}

.form2str <- function (formula) 
{
    if (inherits(formula, "formula")) {
        formula = as.character(formula)
        if (length(formula) == 3) 
            formula = paste(formula[c(2, 1, 3)], collapse = " ")
        else formula = paste(formula, collapse = " ")
    }
    formula
}


code.design <- function(design){
      ## function that applies coding information for subsequent application of rsm and its methods
      if (!"design" %in% class(design)) stop("code.design is applicable to class design objects only.")
      di <- design.info(design)
#! --- RVL suggestion ---
# (coded.data `undesign`s stuff now for fear of conflicts. Guess the opposite happens)
      dn <- attr(design, "desnum")
      ro <- attr(design, "run.order")
      if (is.null(di$coding)) {
          warning("code.design did not change anything,\nbecause the design does not contain coding information.")
          return(design)
          }
      else{
        design.info(design)$coding <- NULL
        hilf <- rep(list(c(-1,1)), di$nfactors)
        names(hilf) <- paste("x",1:di$nfactors, sep="")
        design.info(design)$factor.names <- hilf
#! --- RVL suggestion ---
        des <- coded.data(design, formulas=di$coding)
        }
# Need to add back info...
      class(des) <- c("design", class(des))
      design.info(des) <- di
      attr(des, "desnum") <- dn
      attr(des, "run.order") <- ro
      des
#! --- end RVL suggestion ---
}

decode.design <- function(design){
      ## function that applies coding information for subsequent application of rsm and its methods
      if (!"design" %in% class(design)) stop("decode.design is applicable to class design objects only.")
      if (!"coded.data" %in% class(design)) stop("decode.design is applicable to class coded.data objects only.")

      di <- design.info(design)
      di$coding <- attr(design, "codings")
      fn <- factor.names(design)
      nm = names(fn)
       for (f in di$coding) {
#! --- RVL suggestion ---
            info = .parse.coding(f)   # WAS parse.coding(f)
#! --- end RVL suggestion ---
            cod = info$names[["coded"]]
            org = info$names[["orig"]]
            if (!is.null(fn[[cod]])) {
                fn[[cod]] = info$const[["divisor"]] * fn[[cod]] +
                  info$const[["center"]]
                nm[nm == cod] = org
            }
        }
      names(fn) <- nm
      di$factor.names <- fn
      design.info(design) <- di
      decode.data(design)
}
