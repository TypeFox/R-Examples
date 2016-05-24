# CHNOSZ/util.expression.R
# write descriptions of chemical species, properties, reactions, conditions
# modified from describe(), axis.label()  20120121 jmd

expr.species <- function(species, state="", log="", value=NULL) {
  # make plotting expressions for chemical formulas
  # that include subscripts, superscripts (if charged)
  # and optionally designations of states +/- loga or logf prefix
  if(length(species) > 1) (stop("more than one species"))
  # the counts of elements in the species:
  # here we don't care too much if an "element" is a real element
  # (listed in thermo$element), so we suppress warnings
  elements <- suppressWarnings(makeup(species))
  # where we'll put the expression
  expr <- ""
  # loop over elements
  for(i in 1:length(elements)) {
    if(names(elements)[i] != 'Z') {
      # append the elemental symbol
      expr <- substitute(paste(a, b), list(a=expr, b=names(elements)[i]))
      # recover the coefficient
      if(elements[i]==1) coeff <- "" else coeff <- elements[i]
      # append the coefficient
      # subscripts within subscripts (log) are too small
      if(log != "") expr <- substitute(a*b, list(a=expr, b=coeff))
      else expr <- substitute(a[b], list(a=expr, b=coeff))
    } else {
      # for charged species, don't show "Z" but do show e.g. "+2"
      coeff <- elements[i]
      if(coeff==-1) coeff <- "-"
      else if(coeff==1) coeff <- "+"
      else if(coeff > 0) coeff <- paste("+", as.character(coeff), sep="")
      # append the coefficient (as a superscript if we're not in a log expression)
      if(log != "") expr <- substitute(a*b, list(a=expr, b=coeff))
      else expr <- substitute(a^b, list(a=expr, b=coeff))
    }
  }
  # write a designation of physical state
  # use the state given in log if it's a gas or neutral aqueous species
  if(log %in% c("g", "gas")) state <- "g"
  else if(!"Z" %in% names(elements) & !missing(log)) state <- log
  if(state != "") {
    # subscript it if we're not in a log expression
    if(log != "") expr <- substitute(a*group('(',italic(b),')'),list(a=expr, b=state))
    else expr <- substitute(a[group('(',italic(b),')')],list(a=expr, b=state))
  }
  # write logarithm of activity or fugacity
  if(log != "") {
    if(log %in% c("aq", "cr", "liq", "cr1", "cr2", "cr3", "cr4")) acity <- "a"
    else if(log %in% c("g", "gas")) acity <- "f"
    else stop(paste("'", log, "' is not a recognized state", sep=""))
    logacity <- substitute(log*italic(a), list(a=acity))
    expr <- substitute(a[b], list(a=logacity, b=expr))
    # write a value if given
    if(!is.null(value)) {
      expr <- substitute(a==b, list(a=expr, b=value))
    }
  }
  return(expr)
}

expr.property <- function(property) {
  # a way to make expressions for various properties
  # e.g. expr.property('DG0r') for standard molal Gibbs 
  # energy change of reaction
  propchar <- s2c(property)
  expr <- ""
  # some special cases
  if(property=="logK") return(quote(log~italic(K)))
  # grepl here b/c diagram() uses "loga.equil" and "loga.basis"
  if(grepl("loga", property)) return(quote(log~italic(a)))
  if(property=="alpha") return(quote(alpha))
  if(property=="Eh") return("Eh")
  if(property=="pH") return("pH")
  if(property=="pe") return("pe")
  if(property=="IS") return("IS")
  if(property=="ZC") return(quote(bar(italic(Z))[C]))
  # process each character in the property abbreviation
  prevchar <- character()
  for(i in 1:length(propchar)) {
    if(i > 1) prevchar <- thischar
    thischar <- propchar[i]
    # unless indicated below, uppercase letters are italicized
    # and lowercase letters are italicized and subscripted
    # (includes f for property of formation and r for property of reaction)
    if(thischar %in% LETTERS) thisexpr <- substitute(italic(a), list(a=thischar))
    else if(thischar %in% letters) thisexpr <- substitute(""[italic(a)], list(a=thischar))
    else thisexpr <- substitute(a, list(a=thischar))
    # D for greek Delta
    # A for bold A (affinity)
    # p for subscript italic P (in Cp)
    # 0 for degree sign (but not immediately following a number e.g. 2.303)
    if(thischar=='D') thisexpr <- substitute(Delta)
    if(thischar=='A') thisexpr <- substitute(bold(A))
    if(thischar=='p') thisexpr <- substitute(a[italic(P)], list(a=""))
    if(thischar=='0' & !can.be.numeric(prevchar)) thisexpr <- substitute(degree)
    # put it together
    expr <- substitute(a*b, list(a=expr, b=thisexpr))
  }
  return(expr)
}

expr.units <- function(property, prefix="", per="mol") {
  # make an expression describing units
  # unless we have match below, there will be no units
  # (e.g., logK, pH, pe)
  expr <- ""
  # A, G, H - energy
  if(grepl("A", property)) expr <- substitute(a, list(a=E.units()))
  if(grepl("G", property)) expr <- substitute(a, list(a=E.units()))
  if(grepl("H", property) & !grepl("pH", property)) expr <- substitute(a, list(a=E.units()))
  # Cp, S - energy (per K)
  if(grepl("Cp", property)) expr <- substitute(a~K^-1, list(a=E.units()))
  if(grepl("S", property)) expr <- substitute(a~K^-1, list(a=E.units()))
  # V - volume
  if(grepl("V", property)) expr <- substitute(a^3, list(a="cm"))
  # E - volume (per K)
  if(grepl("E", property)) expr <- substitute(a^3~K^-1, list(a="cm"))
  # P - pressure
  if(grepl("P", property)) expr <- substitute(a, list(a=P.units()))
  # T - temperature
  if(grepl("T", property)) {
    expr <- substitute(a, list(a=T.units()))
    # add a degree sign for Celsius
    if(T.units()=="C") expr <- substitute(degree*a, list(a=expr))
  }
  # Eh - electrical potential
  if(grepl("Eh", property)) expr <- substitute(a, list(a="volt"))
  # IS - ionic strength
  if(grepl("IS", property)) expr <- substitute(a, list(a=mol~kg^-1))
  if(!expr=="") {
    # add prefix if appropriate
    if(!prefix=="") expr <- substitute(a*b, list(a=prefix, b=expr))
    # add mol^-1 if appropriate
    if(!any(sapply(c("P", "T", "Eh", "IS"), function(x) grepl(x, property))))
      expr <- substitute(a~b^-1, list(a=expr, b=per))
  }
  return(expr)
}

axis.label <- function(label, units=NULL, basis=get("thermo")$basis, prefix="") {
  # make a formatted axis label from a generic description
  # it can be a chemical property, condition, or chemical activity in the system
  # if the label matches one of the basis species
  # or if the state is specified, it's a chemical activity
  # 20090826: just return the argument if a comma is already present
  # (it's good for custom labels that shouldn't be italicized)
  if(grepl(",", label)) return(label)
  if(label %in% rownames(basis)) {
    # 20090215: the state this basis species is in
    state <- basis$state[match(label, rownames(basis))]
    # get the formatted label
    desc <- expr.species(label, log=state)
  } else {
    # the label is for a chemical property or condition
    # make the label by putting a comma between the property and the units
    property <- expr.property(label)
    if(is.null(units)) units <- expr.units(label, prefix=prefix)
    # no comma needed if there are no units
    if(units=="") desc <- substitute(a, list(a=property))
    else desc <- substitute(list(a, b), list(a=property, b=units))
  }
  # done!
  return(desc)
}

describe.basis <- function(basis=get("thermo")$basis, ibasis=1:nrow(basis), digits=1, oneline=FALSE) {
  # make expressions for the chemical activities/fugacities of the basis species
  propexpr <- valexpr <- character()
  for(i in ibasis) {
    # propexpr is logarithm of activity or fugacity
    propexpr <- c(propexpr, expr.species(rownames(basis)[i], log=basis$state[i]))
    # we have an as.numeric here in case the basis$logact is character
    # (by inclusion of a buffer for one of the other basis species)
    valexpr <- c(valexpr, format(round(as.numeric(basis$logact[i]), digits), nsmall=digits))
  }
  # write an equals sign between the property and value
  desc <- character()
  for(i in seq_along(propexpr)) {
    thisdesc <- substitute(a==b, list(a=propexpr[[i]], b=valexpr[[i]]))
    if(oneline) {
      # put all the property/value equations on one line, separated by commas
      if(i==1) desc <- substitute(a, list(a=thisdesc))
      else desc <- substitute(list(a, b), list(a=desc, b=thisdesc))
    } else desc <- c(desc, thisdesc)
  }
  return(as.expression(desc))
}

describe.property <- function(property=NULL, value=NULL, digits=1, oneline=FALSE, ret.val=FALSE) {
  # make expressions for pressure, temperature, other conditions
  if(is.null(property) | is.null(value)) stop("property or value is NULL")
  propexpr <- valexpr <- character()
  for(i in 1:length(property)) {
    propexpr <- c(propexpr, expr.property(property[i]))
    thisvalue <- format(round(value[i], digits), nsmall=digits)
    thisunits <- expr.units(property[i])
    thisvalexpr <- substitute(a~b, list(a=thisvalue, b=thisunits))
    valexpr <- c(valexpr, as.expression(thisvalexpr))
  } 
  # with ret.val=TRUE, return just the value with the units (e.g. 55 degC)
  if(ret.val) return(valexpr)
  # write an equals sign between the property and value
  desc <- character()
  for(i in seq_along(propexpr)) {
    thisdesc <- substitute(a==b, list(a=propexpr[[i]], b=valexpr[[i]]))
    if(oneline) {
      # put all the property/value equations on one line, separated by commas
      if(i==1) desc <- substitute(a, list(a=thisdesc))
      else desc <- substitute(list(a, b), list(a=desc, b=thisdesc))
    } else desc <- c(desc, thisdesc)
  }
  return(as.expression(desc))
}

describe.reaction <- function(reaction, iname=numeric(), states=NULL) {
  # make an expression describing the reaction that is
  # the 'reaction' part of subcrt() output
  reactexpr <- prodexpr <- character()
  # loop over the species in the reaction
  for(i in 1:nrow(reaction)) {
    # get the name or the chemical formula of the species
    if(i %in% iname) species <- reaction$name[i]
    else {
      # should the chemical formula have a state?
      if(identical(states,"all")) species <- expr.species(reaction$formula[i], state=reaction$state[i])
      else species <- expr.species(reaction$formula[i])
    }
    # get the absolute value of the reaction coefficient
    abscoeff <- abs(reaction$coeff[i])
    # put the coefficient in if it's not 1
    if(abscoeff==1) coeffspec <- species
    else {
      # we put in some space if the coefficient comes before a name
      if(i %in% iname) coeffspec <- substitute(a~b, list(a=abscoeff, b=species))
      else coeffspec <- substitute(a*b, list(a=abscoeff, b=species))
    }
    # is it a reactant or product?
    if(reaction$coeff[i] < 0) {
      if(length(reactexpr)==0) reactexpr <- substitute(a, list(a=coeffspec))
      else reactexpr <- substitute(a+b, list(a=reactexpr, b=coeffspec))
    } else {
      # it's a product
      if(length(prodexpr)==0) prodexpr <- substitute(a, list(a=coeffspec))
      else prodexpr <- substitute(a+b, list(a=prodexpr, b=coeffspec))
    }
  }
  # put an equals sign between reactants and products
  desc <- substitute(a==b, list(a=reactexpr, b=prodexpr))
  return(desc)
}

