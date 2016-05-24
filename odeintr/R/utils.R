# Copyright Timothy H. Keitt 2015

read_template = function(name)
{
  tn = paste0(name, ".cpp")
  fn = system.file(file.path("templates", tn),
                   package = "odeintr", mustWork = TRUE)
  con = file(fn)
  res = readLines(con)
  close(con)
  paste0(res, collapse = "\n")
}

vectorize_1d_sys = function(sys)
{
  if (grepl("\\[\\s*\\d+\\s*\\]", sys)) return(sys)
  sys = gsub("\\bx\\b", "x[0]", sys)
  sys = gsub("\\bdxdt\\b", "dxdt[0]", sys)
  return(sys)
}

make_stepper_constr = function(method, atol, rtol)
{
  stepper_constr = "stepper_type()"
  if (grepl("euler_|rk4_|rk54_i|rk78_i|bs_|bsd_", method))
    stop("Invalid integration method")
  if (grepl("_a$", method))
    stepper_constr = paste0("odeint::make_controlled(", atol, ", ", rtol, ", stepper_type())")
  if (grepl("_i$", method))
    stepper_constr = paste0("odeint::make_dense_output(", atol, ", ", rtol, ", stepper_type())") 
  return(stepper_constr)
}

make_stepper_type = function(stepper)
{
  stepper = sub("_i$|_a$", "", stepper)
  if (grepl("ab|am|abm", stepper))
  {
    steps = as.integer(sub("ab|am|abm([0-9]+)", "\\1", stepper))
    stepper = sub("(ab|am|abm)[0-9]+", "\\1", stepper)
  }
  switch(stepper,
         euler = "euler<state_type>",
         rk4 = "runge_kutta4<state_type>",
         rk54 = "runge_kutta_cash_karp54<state_type>",
         rk5 = "runge_kutta_dopri5<state_type>",
         rk78 = "runge_kutta_fehlberg78<state_type>",
         ab = paste0("adams_bashforth<", steps, ", state_type>"),
         am = paste0("adams_moulton<", steps, ", state_type>"),
         abm = paste0("adams_bashforth_moulton<", steps, ", state_type>"),
         bs = "bulirsch_stoer<state_type>",
         bsd = "bulirsch_stoer_dense_out<state_type>",
         paste0(stepper, "<state_type>"))
}

get_sys_dim = function(x)
{
  x = gsub("\\[\\s*(\\d+)\\s*\\]", "\\[\\1\\]", x)
  matches = gregexpr("dxdt\\[\\d+\\]", x)
  lens = attr(matches[[1]], "match.length") - 7L
  starts = unlist(matches) + 5L
  indices = rep(NA, length(starts))
  for (i in seq(along = indices))
    indices[i] = substr(x, starts[i], starts[i] + lens)
  return(max(as.integer(indices)) + 1L)
}

disable_asserts = function(makevars)
{
  con = pipe(paste("R CMD config CPPFLAGS"))
  flags = readLines(con); close(con)
  if (!is.finite(pmatch("-DNDEBUG", flags)))
    flags = paste(flags, "-DNDEBUG")
  flags = gsub("^\\s+|\\s+$", "", flags)
  cat(paste0("CPPFLAGS=", flags, "\n"),
      file = makevars, append = TRUE)
}

substitute_opt_level = function(flags, level, omit.debug)
{
  flags = gsub("-O\\d+", paste0("-O", level), flags)
  if (omit.debug) flags = gsub("\\s*-g\\s*", "", flags)
  flags = gsub("^\\s+|\\s+$", "", flags)
  return(flags)
}

process_flags = function(name, level, omit.debug)
{
  con = pipe(paste("R CMD config", name))
  flags = readLines(con); close(con)
  flags = substitute_opt_level(flags, level, omit.debug)
  paste0(name, "=", flags)
}

Jacobian1 = function(f)
{
  sep = ".."
  vn = names(formals(f))[1]
  e = body(f)
  for (i in length(e))
  {
    es = deparse(e[[i]])
    if (grepl(paste0("\\b", vn, "\\b"), es))
    {
      es = gsub("\\[\\s*(\\d+)\\s*\\]", paste0(sep, "\\1"), es)
      de = deparse(D(parse(text = es), vn))
      de = gsub(paste0("\\bx", sep, "(\\d+)"), "x\\[\\1\\]", de)
      e[[i]] = parse(text = de)
    }
  }
  res = function() NULL
  formals(res) = formals(f)
  body(res) = e
  return(res)
}

Jacobian2 = function(code, sys_dim = -1)
{
  sep = ".."
  code = paste0(code, collapse = ";")
  if (sys_dim < 1)
    sys_dim = get_sys_dim(code)
  if (is.na(sys_dim))
  {
    code = vectorize_1d_sys(code)
    sys_dim = 1L
  }
  code = gsub("^\\s+|\\s+$", "", unlist(strsplit(code, ";")))
  code = code[nzchar(code) != 0]
  i = unlist(lapply(code, function(x)
    as.numeric(sub("\\bdxdt\\[\\s*(\\d+)\\s*\\].*", "\\1", x))))
  code = code[order(i)]
  g = function(j, i, rhs)
  {
    var = paste0("x", sep, j - 1)
    deriv = D(parse(text = rhs), var)
    deriv = deparse(deriv)
    deriv = gsub(paste0("\\bx", sep, "(\\d+)"), "x\\[\\1\\]", deriv)
    deriv = paste0("J(", i - 1, ", ", j - 1, ") = ", deriv, ";") 
  }
  f = function(i)
  {
    rhs = sub("\\bdxdt\\[\\s*\\d+\\s*\\]\\s*=\\s*(.*)", "\\1", code[i])
    rhs = gsub("\\[\\s*(\\d+)\\s*\\]", paste0(sep, "\\1"), rhs)
    return(unlist(lapply(1:sys_dim, g, i, rhs)))
  }
  paste0(unlist(lapply(1:sys_dim, f)), collapse = "\n")
}

handle_pars = function(pars, const = FALSE)
{
  switch(mode(pars),
         numeric =
         {
           pn = names(pars)
           if (!is.null(pn) && !all(nzchar(pn)))
             stop("All parameters must be named")
           if (is.null(pn))
             make_vec_pars_decl(pars[[1]])
           else
             if (const)
              list(globals = make_init_pars_decl(pars, TRUE))
             else
               list(globals = make_init_pars_decl(pars),
                    setter = make_init_pars_setter(pars),
                    getter = make_init_pars_getter(pars))
         },
         character =
           list(globals = make_pars_decl(pars),
                setter = make_pars_setter(pars),
                getter = make_pars_getter(pars)),
         stop("Invalid parameter specification"))
}

make_init_pars_decl = function(pars, const = FALSE)
{
  res = paste(names(pars), "=", pars, collapse = ", ")
  res = paste0("double ", res, ";")
  if (const) res = paste("const", res)
  res
}

make_init_pars_setter = function(pars)
{
  pn = names(pars)
  body = paste0("odeintr::", pn, " = ", pn, collapse = ";\n")
  pars = paste0("double ", pn, " = ", pars, collapse = ", ")
  code = read_template("pars_setter_template")
  code = sub("__PARS__", pars, code)
  code = sub("__BODY__", body, code)
  return(paste0(code, collapse = "\n"))
}

make_init_pars_getter = function(pars)
{
  pn = names(pars)
  body = paste0("out[\"", pn, "\"] = odeintr::", pn, collapse = ";\n")
  code = read_template("pars_getter_template")
  code = sub("__BODY__", body, code)
  return(paste0(code, collapse = "\n"))
}

make_vec_pars_decl = function(N)
{
  setter = read_template("pars_vec_setter_template")
  getter = read_template("pars_vec_getter_template")
  list(globals = paste0("std::array<double, ", N, "> pars;"),
       setter = paste0(setter, collapse = "\n"),
       getter = paste0(getter, collapse = "\n"))  
}

make_pars_decl = function(pars)
{
  res = paste(pars, collapse = ", ")
  return(paste0("double ", res, ";"))
}

make_pars_setter = function(pars)
{
  body = paste0("odeintr::", pars, " = ", pars, collapse = ";\n")
  pars = paste0("double ", pars, collapse = ", ")
  code = read_template("pars_setter_template")
  code = sub("__PARS__", pars, code)
  code = sub("__BODY__", body, code)
  return(paste0(code, collapse = "\n"))
}

make_pars_getter = function(pars)
{
  body = paste0("out[\"", pars, "\"] = odeintr::", pars, collapse = ";\n")
  code = read_template("pars_getter_template")
  code = sub("__BODY__", body, code)
  return(paste0(code, collapse = "\n"))
}

proc_output = function(res)
{
  if (length(res[[1]]) == 0) return(NULL)
  x = res[[2]]; out = list(res[[1]])
  if (any(diff(sapply(x, length)) != 0)) return(res)
  out = append(out, rw2cw(x))
  xnames = names(x[[1]])
  if (is.null(xnames) || length(xnames) != length(x[[1]]))
    xnames = paste0("X", 1:length(x[[1]]))
  names(out) = c("Time", xnames)
  attr(out, "row.names") = c(NA, -length(out[[1]]))
  class(out) = "data.frame"
  return(out)
}

rw2cw = function(x)
{
  lapply(lapply(1:length(x[[1]]),
                function(i) lapply(x, function(a) a[[i]])),
         unlist)
}

