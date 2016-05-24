# how are the overrides structured?
  # list(
  #   GtkWidget = list(method1 = function() {...},
  #                    method2 = function() {...},
  #                    ...
  #               ),
  #   GtkContainer = list(...),
  #   ...
  # )
  #

.reserved <- c(".props", ".prop_overrides", ".initialize", ".signals",
    ".public", ".protected", ".private")

# FIXME: Should we really require that interfaces be specified only within
# the class definition - no explicit argument?
gClass <- function(name, parent = "GObject", ..., abstract = FALSE)
{
  virtuals <- as.list(.virtuals)
  
  # make sure name is valid
  
  name <- as.character(name)
  if (!length(name) || nchar(name) < 3)
    stop("Type name must be coerceable to a string and be at least 3 characters long")
  invalid_chars <- sub("^[a-zA-Z_][a-zA-Z+\\-_0-9]*", "", name)
  if (nchar(invalid_chars))
    stop("The characters \"", invalid_chars, "\" are invalid. Names must start with",
    " a letter or '_' and contain only letters, numbers, or any of '+-_'")
  
  # validate parent type
  
  parent <- as.character(parent)
  full_hierarchy <- gTypeGetAncestors(parent)
  if (!length(full_hierarchy) || full_hierarchy[length(full_hierarchy)] != "GObject")
    stop("Parent type ", parent, " is not a GObject-derived type")
  
  # process class definition

  class_list <- list(...)
  missing_method_lists <- setdiff(c(".public", ".protected", ".private"),
                                  names(class_list))
  class_list[missing_method_lists] <-
    list(structure(list(), names = character()))
  types <- class_list[!(names(class_list) %in% .reserved)]
  known_types <- names(types) %in% names(virtuals)
  if (any(!known_types))
    warning("The types ", paste(names(types)[!known_types], collapse=","), 
      " are not recognized")
  types <- types[known_types]
  ancestors <- lapply(names(types), gTypeGetAncestors)
  interfaces <- as.logical(sapply(ancestors, function(hierarchy) 
    hierarchy[length(hierarchy)] == "GInterface"))
  classes <- types[!interfaces]
  interface_names <- names(types)[interfaces]
  
  # try to detect problems in class definition
  
  ancestors <- ancestors[!interfaces]
  ancestor_lengths <- sapply(ancestors, length)
  parallel_hierarchies <- duplicated(ancestor_lengths)
  if (any(parallel_hierarchies))
    stop("The following classes have parallel hierarchies: ", 
      paste(sapply(ancestor_lengths[parallel_hierarchies], function(dup) 
        paste("(", paste(classes[ancestor_lengths == dup], collapse=","), ")", sep=""),
      collapse=",")))
  
  inconsistents <- !(unlist(ancestors) %in% full_hierarchy)
  if (any(inconsistents))
    stop("Inconsistent hierarchy detected. The classes ", 
      paste(unique(unlist(ancestors)[inconsistents]), collapse=","),
      " are not present in the full hierarchy: ", paste(full_hierarchy, collapse=","))
  
  # fill in the gaps
  
  missing_classes <- which(!(full_hierarchy %in% names(classes)))
  types[full_hierarchy[missing_classes]] <- replicate(length(missing_classes), list())
  
  # don't want to fill in the virtuals for SGObject parents
  sapply(names(types), function(type_name)
    if (!("SGObject" %in% gTypeGetInterfaces(type_name)))
      types[[type_name]] <<- types[[type_name]][virtuals[[type_name]]])
  
  # check init function
  
  init <- class_list$.init
  if (!is.null(init))
    init <- as.function(init)
  
  # check property paramspecs
  
  props <- lapply(class_list$.props, as.GParamSpec)
  
  # and get any overriden properties
  
  prop_overrides <- as.character(class_list$.prop_overrides)
  
  # check signals
  
  signal_fields <- c("name", "param_types", "ret_type", "flags")
  signals <- lapply(class_list$.signals, function(signal) {
    signal <- as.struct(signal, "GSignal", signal_fields)
    name <- as.character(signal[[1]])
    invalid_chars <- sub("^[a-zA-Z][a-zA-Z_-]*", "", name)
    if (nchar(invalid_chars))
      stop("Invalid signal name: ", name, ". Signal names must start with a letter",
        " and contain only letters, '-', or '_'")
    signal[[1]] <- name
    flags <- signal[[4]]
    if (is.null(flags))
      flags <- c("run-last", "action")
    if (is.character(flags))
      flags <- sum(GSignalFlags[flags])
    signal[[4]] <- as.numeric(flags)
    if (is.null(signal[[3]]))
      signal[[3]] <- "void"
    signal[[3]] <- as.GType(signal[[3]])
    signal[[2]] <- lapply(signal[[2]], as.GType)
    signal
  })
  
  whichFuncs <- function(env, syms) 
    as.logical(unlist(sapply(mget(syms, env), is.function)))
  
  # get the public members (env locked, fields locked, functions locked)
  # public env is static, so always inherits from parent
  public_parent <- emptyenv()
  if ("SGObject" %in% gTypeGetInterfaces(parent))
    public_parent <- .Call("S_g_object_type_get_public_env", parent, PACKAGE = "RGtk2")
  publics <- list2env(class_list$.public, parent = public_parent)
  public_syms <- ls(publics)
  public_which_funcs <- whichFuncs(publics, public_syms)
  # this can be retrieved by R, so lock it now
  lockEnvironment(publics, TRUE)
  
  # get the protected members (env locked, fields unlocked, functions locked)
  # protected is attached to parent after being cloned
  protecteds <- list2env(class_list$.protected)
  protected_syms <- ls(protecteds)
  protected_which_funcs <- whichFuncs(protecteds, protected_syms)
  
  # get the private members (env unlocked, fields and methods unlocked)
  # private env inherits from protected after cloning
  privates <- list2env(class_list$.private)
  private_syms <- ls(privates)
  
  # check for conflicts between declared members
  
  syms <- c(private_syms, protected_syms, public_syms)
  sym_dups <- duplicated(syms)
  if (any(sym_dups))
    stop("Duplicate symbols in class definition: ", 
      paste(unique(syms[sym_dups]), collapse = ", "))
  
  # check for conflicts with ancestors
  
  # just returns NULL when it can't find something
  mget_null <- function(x, envir) 
    mget(x, envir, ifnotfound = vector("list", length(x)))
    
  ancestor_syms <- unlist(c(mget_null(names(types), .virtuals), 
    mget_null(names(types), .fields)))
  ancestor_conflicts <- ancestor_syms %in% syms
  if (any(ancestor_conflicts))
    stop("Declared symbols already declared in parent: ", 
      paste(ancestor_syms[ancestor_conflicts], collapse = ", "))
  
  # assume all public and protected methods are "virtuals"
  my_virtuals <- list(c(public_syms[public_which_funcs], 
    protected_syms[protected_which_funcs]))
  names(my_virtuals) <- name
  registerVirtuals(my_virtuals)
  
  # remember all public and protected fields to avoid conflicts
  my_fields <- list(c(public_syms[!public_which_funcs],
    protected_syms[!protected_which_funcs]))
  names(my_fields) <- name
  .registerFields(my_fields)
  
  # create new type that extends specified class and implements specified interfaces
  
  get_class_init_funcs <- function(class_name) paste("S", sapply(class_name, .collapseClassName),
    "class_init", sep="_")
  class_init_funcs <- get_class_init_funcs(full_hierarchy)
  parent_class_init <- class_init_funcs[sapply(class_init_funcs, is.loaded)][1]
  class_init_sym <- getNativeSymbolInfo(parent_class_init)$address
  interface_init_syms <- NULL
  if (length(interface_names))
    interface_init_syms <- sapply(get_class_init_funcs(interface_names),
      function(interface_name) getNativeSymbolInfo(interface_name)$address)
  
  class_env <- list2env(types)
  assign(".initialize", init, class_env)
  
  # add public, protected, and private to the class env
  # they need to be cloned and given the corresponding cloned environments
  # from the parent class (except private) during instantiation
  assign(".public", publics, class_env)
  assign(".protected", protecteds, class_env)
  assign(".private", privates, class_env)
  
  # put properties in class_env too, so that they can be registered during
  # class initialization
  assign(".props", props, class_env)
  assign(".prop_overrides", prop_overrides, class_env)
  
  .RGtkCall("S_gobject_class_new", name, parent, interface_names, 
    class_init_sym, interface_init_syms, class_env, signals, abstract)
}

registerVirtuals <- function(virtuals)
{
  list2env(virtuals, .virtuals)
}

.registerFields <- function(fields)
{
  list2env(fields, .fields)
}

gObjectParentClass <- function(obj)
{
  gTypeGetClass(class(obj)[2])
}

assignProp <- function(obj, pspec, value)
{
  stopifnot(implements(obj, "SGObject"))
  assign(pspec$name, value, attr(obj, ".private"))
}
getProp <- function(obj, pspec)
{
  stopifnot(implements(obj, "SGObject"))
  get(pspec$name, attr(obj, ".private"))
}

# takes vector of class env's starting from root, clones the protected env's,
# establishes the inheritance structure, and then clones the private env
# for the leaf type, inheriting from its protected env
# finally, needs to copy overrides from static env to public/protected envs
# the functions are static, but this way we just worry about 2 environments in R
.instanceEnv <- function(hierarchy)
{
  static <- hierarchy[[length(hierarchy)]]
  prot <- emptyenv()
  sapply(hierarchy, function(env) {
    prot <<- list2env(as.list(env$.protected), parent = prot)
  })
  priv <- list2env(as.list(static$.private), parent = prot)
  pub <- static$.public
  assign(".public", pub, static) # ????
  # if something in static exists in public or protected (parents), copy it over
  static_syms <- ls(static)
  static_syms <- setdiff(static_syms, .reserved)
  list2env(mget(intersect(static_syms, ls(pub)), static), pub)
  list2env(mget(intersect(static_syms, ls(priv)), static), priv)
  # note that C overrides are not copied here, and they shouldn't be
  # need to lock up the envs now (public completely, protected just funcs)
  lockEnvironment(pub, TRUE)
  lockEnvironment(prot)
  prot_syms <- ls(prot)
  prot_funcs <- as.logical(sapply(mget(prot_syms, prot), is.function))
  sapply(prot_syms[prot_funcs], lockBinding, prot)
  priv
}

