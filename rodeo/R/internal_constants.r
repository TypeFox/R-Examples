
rodeoConst= list(

  # Regular expressions to detect valid identifiers
  identifierPatterns= list(
    core="[a-zA-Z]([a-zA-Z0-9_])*",
    before= "(^|[^a-zA-Z0-9_])",
    after= "([^a-zA-Z0-9_]|$)"
  ),

  # Identifiers used in generated code
  genIdent= list(
    # vectors with actual values of vars, pars, rates, derivatives
    vec=     c(vars="var",    pars="par",    pros="pro", drvs="dydt"),
    # index vectors for the 0D case (or 1st level, respectively)
    ivec0D=  c(vars="ivar0D", pars="ipar0D", pros="ipro0D"),
    # index vectors for the current level
    ivec=    c(vars="ivar",   pars="ipar",   pros="ipro"),
    # (pointers to) indices of variables in first level
    ilist0D= c(vars="v0D",    pars="p0D",    pros="r0D"),
    # (pointers to) indices of variables in current level
    ilist=   c(vars="v",      pars="p",      pros="r"),
    # dimensions
    len=     c(vars="NVAR",   pars="NPAR",   pros="NPRO"),
    # number of levels and index of current level
    nLevels="NLVL", iLevel= "level"
  ),

  reservedNames= c(
    time="time",     # external time
    left="left",     # reference to left element in a 1D model
    right="right"    # reference to right element in a 1D model
  ),

  lang= c(r="r", fortran="f95")

) # End of list

