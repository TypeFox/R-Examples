###################################
### generic registry infrastructure

### IDEA: use lexical scope with nested functions to create an
### S3-"object" that exposes the data structure only through accessor functions.

.FUNCall <- function(f) function(...) f(...)

registry <-
function(index_field = "names", entry_class = NULL,
         validity_FUN = NULL, registry_class = NULL,
         ignore_case = TRUE)
{
### ATTRIBUTES
    ## repository
    DATA <- META <- list()

    ## permissions
    PERMISSIONS <- c(set_entries = TRUE, modify_entries = TRUE,
                     delete_entries = TRUE, set_fields = TRUE)
    SEALED_FIELDS <- SEALED_ENTRIES <- character(0)


### METHODS (PRIVATE)
    ## helper functions
    .field_exists <-
    function(name)
        name %in% .get_field_names()

    .make_field <-
    function(default = NA, type = NA, is_mandatory = FALSE,
             is_modifiable = TRUE, validity_FUN = NULL)
        structure(list(type = type,
                       default = default,
                       is_mandatory = is_mandatory,
                       is_modifiable = is_modifiable,
                       validity_FUN = validity_FUN),
                  class = "registry_field")

    .make_entry <- function(l)
    {
        ## sort
        l <- l[c(index_field, setdiff(.get_field_names(), index_field))]

        ## return object (possibly inheriting from entry_class)
        structure(l, class = c(entry_class, "registry_entry"))
    }

    .get_mandatory_fields <-
    function()
        names(META)[sapply(META, function(i) i$is_mandatory)]

    .get_field_defaults <-
    function()
        lapply(META, function(i) i$default)

    .get_entry_index <-
    function(name, stop_if_missing = TRUE) {
        ## returns the index of the first exact match (modulo case):
        index <- if (ignore_case)
            sapply(DATA,
                   function(i) toupper(name) %in% toupper(i[[index_field]]))
        else
            sapply(DATA,
                   function(i) name %in% i[[index_field]])
        if (!any(index)) {
            if (stop_if_missing)
                stop(paste("Entry", dQuote(name), "not in registry."))
            else
                return(NULL)
        }
        which(index)[1]
    }

    .check_value <-
    function(field_name, field, value)
    {
        ## Note we do not check NA entries because this may by set automatically
        if (!is.function(value) && !is.na(value[1])) {
            ## check class / list of alternatives, if any
            if (!any(is.na(field$type))) {
                ## check list of alternatives
                if (length(field$type) > 1) {
                    if (!is.character(value) || !value %in% field$type)
                        stop(paste("Possible values for", dQuote(field_name), "are:",
                                   paste(field$type, collapse = ", ")))
                ## check class
                } else if (!inherits(value, field$type)) {
                    stop(paste("Field", dQuote(field_name), "does not inherit from class", field$type))
                }
            }

            ## apply validity function, if any
            if (!is.null(field$validity_FUN))
                do.call(field$validity_FUN, list(value))
        }
    }

    .check_for_unknown_fields <-
    function(n)
    {
        ## check for fields not in repository
        missing <- !.field_exists(n)
        if (any(missing))
            stop(paste("Field(s) not in repository:", paste(n[missing], collapse = ", ")))
    }

### METHODS (PUBLIC)
    ## field accessors
    .entry_exists <-
    function(name)
        if (ignore_case)
            toupper(name) %in% toupper(sapply(DATA, function(i) i[[index_field]]))
        else
            name %in% sapply(DATA, function(i) i[[index_field]])

    .get_field <-
    function(name)
    {
        if (!.field_exists(name))
            stop(paste("Field", dQuote(name), "not in registry."))

        META[[name]]
    }

    .get_fields <-
    function()
        META

    .get_field_names <-
    function()
        names(META)

    .set_field <-
    function(name, default = NA, type = NA, is_mandatory = FALSE,
             is_modifiable = TRUE, validity_FUN = NULL)
    {
        ## check permissions
        if (!PERMISSIONS["set_fields"])
            stop("Setting of fields not allowed.")

        ## check for double entries
        if (.field_exists(name))
            stop(paste("Field", dQuote(name), "already in registry."))

        ## possibly, infer type from argment
        if (!is.na(type) && !(is.character(type)))
            type <- class(type)

        ## check mandatory fields
        if (is_mandatory && !is.na(default))
            stop("Mandatory fields should have no default.")

        ## create field entry
        field <- .make_field(type = type,
                             default = default,
                             is_mandatory = is_mandatory,
                             validity_FUN = validity_FUN)

        ## check validity of default
        .check_value("default", field, default)

        ## add field to meta data
        META <<- c(META, list(field))
        names(META)[length(META)] <<- name

        ## add (missing) fields to data entries
        DATA <<- lapply(DATA, function(i) {i[[name]] <- default; i})
    }

    .n_of_entries <-
    function()
        length(DATA)

    ## entry accessors
    .set_entry <-
    function(...)
    {
        ## check permissions
        if (!PERMISSIONS["set_entries"])
            stop("Setting of entries not allowed.")

        ## parameter handling
        l <- list(...)
        n <- names(l)

        .check_for_unknown_fields(n)

        ## check for mandatory fields
        mandatory_fields <- .get_mandatory_fields()
        missing_mandatory_fields <- !mandatory_fields %in% n
        if (any(missing_mandatory_fields))
            stop(paste("The following fields are mandatory, but missing:",
                       paste(mandatory_fields[missing_mandatory_fields], collapse = ", ")))

        ## check for double entries
        for (i in l[[index_field]])
            if (.entry_exists(i))
                stop(paste("Entry", dQuote(i), "already in registry."))

        ## check defaults and set values, if needed
        field_defaults    <- .get_field_defaults()
        default_fields    <- names(field_defaults)
        missing_fields    <- setdiff(default_fields, n)
        l[missing_fields] <- field_defaults[missing_fields]

        ## check field types, and apply field check function, if any.
        for (f in n) {
            meta <- .get_field(f)
            .check_value(f, .get_field(f), l[[f]])
        }

        ## apply entry check function
        if (!is.null(validity_FUN))
            do.call(validity_FUN, list(l))

        ## add entry
        entry <- .make_entry(l)
        DATA <<- c(DATA, list(entry))
        names(DATA)[length(DATA)] <<- l[[index_field]][1]
    }

    .get_entries <-
    function(names = NULL, pattern = NULL) {
        ## fix search
        if (!is.null(names)) {
            if (ignore_case)
                DATA[intersect(toupper(names), toupper(names(DATA)))]
            else
                DATA[intersect(names, names(DATA))]
        ## grep search
        } else if (!is.null(pattern)) {
            pattern_in_entry <-
                function(x) any(sapply(x, function(i) is.character(i)
                                       && length(grep(pattern, i) > 0)))
            DATA[sapply(DATA, pattern_in_entry)]
        ## else: return all entries
        } else DATA

    }

    .get_entry_names <-
    function()
    {
        if (length(DATA) < 1)
            character(0)
        else
            names(DATA)
    }

    .get_entry <-
    function(name, stop_if_missing = TRUE)
    {
        index <- .get_entry_index(name, stop_if_missing)
        if (is.null(index))
            return(NULL)
        DATA[[index]]
    }

    .delete_entry <-
    function(name)
    {
        ## check permissions
        if (!PERMISSIONS["delete_entries"])
            stop("Deletion of entries not allowed.")

        ## fetch entry index (this also checks if the entry exists)
        entry_index <- .get_entry_index(name)

        ## check sealed entries
        if (name %in% SEALED_ENTRIES)
            stop(paste("Deletion of entry", dQuote(name), "not allowed."))

        ## delete it
        DATA[entry_index] <<- NULL
    }


    .modify_entry <-
    function(...)
    {
        ## check permissions
        if (!PERMISSIONS["modify_entries"])
            stop("Modifying of entries not allowed.")

        ## parameter handling
        l <- list(...)
        n <- names(l)

        ## check for index field
        if (!index_field %in% n)
            stop(paste("Index field", dQuote(index_field), "missing."))

        .check_for_unknown_fields(n)

        ## determine entry name
        name <- l[[index_field]][1]

        ## fetch entry index (this also checks if the entry exists)
        entry_index <- .get_entry_index(name)

        ## fetch entry
        entry <- DATA[[entry_index]]
        name <- entry[[index_field]][1]

        for (field in setdiff(n, index_field)) {
            ## check if field is modifiable
            field_entry <- .get_field(field)
            if (!field_entry$is_modifiable)
                stop(paste("Field", dQuote(field), "is not modifiable."))

            ## check if entry and field are sealed
            if ((name %in% SEALED_ENTRIES) && (field %in% SEALED_FIELDS))
                stop(paste("Modification of field", dQuote(field),
                           "in entry", dQuote(name), "not allowed."))

            ## check new value
            value <- l[[field]]
            .check_value(field, field_entry, value)

            ## modify entry locally
            entry[[field]] <- value
        }

        ## apply entry check function
        if (!is.null(validity_FUN))
            do.call(validity_FUN, list(entry))

        ## modify entry in registry
        DATA[entry_index] <<- list(entry)
    }

    ## get all entries for one field
    .get_field_entries <-
    function(field, unlist = TRUE)
    {
        if (!.field_exists(field))
            stop(paste("Field", dQuote(field), "not in registry."))

        ret <- lapply(DATA, function(i) i[[field]])
        if (unlist)
            unlist(ret)
        else
            ret
    }

    ## permission getters/setters
    .get_permissions <-
    function()
        PERMISSIONS

    .restrict_permissions <-
    function(set_entries = TRUE, modify_entries = TRUE,
             delete_entries = TRUE, set_fields = TRUE)
    {
        PERMISSIONS["set_entries"] <<- PERMISSIONS["set_entries"] && set_entries
        PERMISSIONS["modify_entries"] <<- PERMISSIONS["modify_entries"] && modify_entries
        PERMISSIONS["delete_entries"] <<- PERMISSIONS["delete_entries"] && delete_entries
        PERMISSIONS["set_fields"] <<- PERMISSIONS["set_fields"] && set_fields
    }

    .seal_entries <-
    function()
    {
        SEALED_ENTRIES <<- .get_entry_names()
        SEALED_FIELDS <<- .get_field_names()
    }

    .get_sealed_entry_names<-
    function()
        SEALED_ENTRIES

    .get_sealed_field_names <-
    function()
        SEALED_FIELDS

### CONSTRUCTOR

    ## create index field
    .set_field(name = index_field, type = "character",
               is_mandatory = TRUE, is_modifiable = FALSE)

    ## return class
    structure(list(get_field = .FUNCall(.get_field),
                   get_fields = .FUNCall(.get_fields),
                   get_field_names = .FUNCall(.get_field_names),
                   set_field = .FUNCall(.set_field),

                   entry_exists = .FUNCall(.entry_exists),
                   get_entry = .FUNCall(.get_entry),
                   get_entries = .FUNCall(.get_entries),
                   get_entry_names = .FUNCall(.get_entry_names),
                   set_entry = .FUNCall(.set_entry),
                   modify_entry = .FUNCall(.modify_entry),
                   delete_entry = .FUNCall(.delete_entry),
                   n_of_entries = .FUNCall(.n_of_entries),

                   get_field_entries = .FUNCall(.get_field_entries),

                   get_permissions = .FUNCall(.get_permissions),
                   restrict_permissions = .FUNCall(.restrict_permissions),
                   seal_entries = .FUNCall(.seal_entries),
                   get_sealed_entry_names = .FUNCall(.get_sealed_entry_names),
                   get_sealed_field_names = .FUNCall(.get_sealed_field_names)
                   ),
              class = c(registry_class, "registry"))
}

"[[.registry" <-
function(x, i)
  x$get_entry(i)

length.registry <-
function(x)
    x$n_of_entries()

print.registry <-
function(x, ...)
{
    l <- x$n_of_entries()
    if (l < 1)
        writeLines(paste("An object of class", dQuote("registry"), "with no entry."))
    else if (l == 1)
        writeLines(paste("An object of class", dQuote("registry"), "with one entry."))
    else
        writeLines(paste("An object of class", dQuote("registry"), "with", l, "entries."))
}

print.registry_field <-
function(x, ...)
    writeLines(formatUL(x, label = names(x), ...))

print.registry_entry <-
function(x, ...)
{
    x <- .functions_to_characters(x)
    x[[1]] <- paste(x[[1]], collapse = ", ")
    writeLines(formatUL(x, label = names(x)))
}

summary.registry <-
function(object, ...)
    as.data.frame(object, ...)

as.data.frame.registry <-
function(x, ...)
    do.call(rbind,
            lapply(x$get_entries(),
                   function(entry) {
                       entry <- .functions_to_characters(entry)
                       data.frame(unclass(entry[-1]), ...)
                   }
                   )
            )

.functions_to_characters <-
function(x)
{
    ## transform function entries into character strings
    funs <- sapply(x, inherits, "function")
    for (field in names(x)[funs])
        x[[field]] <- paste(format(x[[field]]), collapse = "")
    x
}
