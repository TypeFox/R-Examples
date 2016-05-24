###################################
### generic registry infrastructure

### IDEA: use lexical scope with nested functions to create an
### S3-"object" that exposes the data structure only through accessor functions.

## creating function
registry <-
    function(entry_class = NULL, registry_class = NULL,
             validity_FUN = NULL, stop_if_missing = FALSE)
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
    function(type = NA, alternatives = NA, default = NA,
             is_mandatory = FALSE, is_modifiable = TRUE, is_key = FALSE,
             validity_FUN = NULL, index_FUN = match_ignorecase,
             index_FUN_args = NULL
             )
        structure(list(type = type,
                       alternatives = alternatives,
                       default = default,
                       is_mandatory = is_mandatory,
                       is_modifiable = is_modifiable,
                       is_key = is_key,
                       validity_FUN = validity_FUN,
                       index_FUN = index_FUN,
                       index_FUN_args = index_FUN_args
                       ),
                  class = "registry_field")


    .make_entry <- function(l)
    {
        ## get index fields
        index_fields <- .get_index_fields()

        ## sort
        l <- l[c(index_fields, setdiff(.get_field_names(), index_fields))]

        ## return object (possibly inheriting from entry_class)
        structure(l, class = c(entry_class, "registry_entry"))
    }

    .get_index_fields <-
    function()
        names(META)[sapply(META, function(i) i$is_key)]

    .get_mandatory_fields <-
    function()
        names(META)[sapply(META, function(i) i$is_mandatory)]

    .get_field_defaults <-
    function()
        lapply(META, function(i) i$default)

    .get_entry_indices <-
    function(key, stop_missing = stop_if_missing, FUN = NULL, ARGS = NULL) {
        ## get index fields
        index_fields <- .get_index_fields()
        l <- length(key)
        n <- names(key)
        if (l != length(index_fields))
            index_fields <- if (is.null(n))
                index_fields[1:l]
            else
                intersect(n, index_fields)

        ## returns the index of the first match
        index <- seq_along(DATA)
        for (index_field in seq_along(index_fields)) {
            if (length(index) < 1)
                break

            if (is.null(FUN)) {
                FUN <- .get_field(index_fields[index_field])$index_FUN
                ARGS <- .get_field(index_fields[index_field])$index_FUN_args
            }
            index <-
                index[sapply(DATA[index],
                        function(i) do.call(FUN,
                                       c(list(key[index_field],
                                              i[[ index_fields[index_field] ]]),
                                         ARGS))
                             )]
        }

        ## nothing found
        if (length(index) < 1) {
            if (stop_missing)
                stop("Entry not in registry.", call. = FALSE)
            else
                return(NULL)
        }

        ## return
        index
    }

    .get_first_entry_index <-
    function(key, stop_missing = stop_if_missing)
        .get_entry_indices(key, stop_missing)[1]

    .check_value <-
    function(field_name, field, value)
    {
        ## Note we do not check NA entries because this may by set automatically
        if (is.object(value) || !is.function(value) && !is.na(value[1])) {
            ## check class / list of alternatives, if any
            if (!any(is.na(field$type))) {
                ## check class
                if (!inherits(value, field$type)) {
                    stop(paste("Field", dQuote(field_name),
                               "does not inherit from:",
                               paste(field$type, collapse = ", ")),
                         call. = FALSE)
                }

                ## check list of alternatives
                if (!any(is.na(field$alternatives))) {
                    if(!any(unlist(lapply(field$alternatives,
                                          identical, value))))
                        stop(paste("Possible values for", dQuote(field_name),
                                   "are:", paste(field$alternatives,
                                                 collapse = ", ")),
                             call. = FALSE)
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
            stop(paste("Field(s) not in repository:",
                       paste(n[missing], collapse = ", ")), call. = FALSE)
    }

### METHODS (PUBLIC)
    ## field accessors
    .get_field <-
    function(name)
    {
        if (!.field_exists(name))
            stop(paste("Field", dQuote(name), "not in registry."),
                 call. = FALSE)

        META[[name]]
    }

    .get_fields <-
    function()
        META

    .get_field_names <-
    function()
        names(META)

    .set_field <-
    function(name,
             type = NA, alternatives = NA, default = NA,
             is_mandatory = NULL, is_modifiable = TRUE, is_key = FALSE,
             validity_FUN = NULL, index_FUN = match_ignorecase, ...
             )
    {
        ## check permissions
        if (!PERMISSIONS["set_fields"])
            stop("Setting of fields not allowed.", call. = FALSE)

        ## check for double entries
        if (.field_exists(name))
            stop(paste("Field", dQuote(name), "already in registry."),
                 call. = FALSE)

        ## possibly, infer type from argment
        if (!is.na(type) && !(is.character(type)))
            type <- class(type)

        ## check mandatory fields
        if (is.null(is_mandatory))
            is_mandatory <- is_key
        if (is_mandatory && !is.na(default))
            stop("Mandatory fields should have no default.", call. = FALSE)
        if (is_key && !is_mandatory)
            stop("Key fields must be mandatory.", call. = FALSE)

        ## create field entry
        field <- .make_field(type = type,
                             alternatives = alternatives,
                             default = default,
                             is_mandatory = is_mandatory,
                             is_modifiable = is_modifiable,
                             is_key = is_key,
                             validity_FUN = validity_FUN,
                             index_FUN = index_FUN,
                             index_FUN_args = list(...)
                             )

        ## check validity of default
        .check_value("default", field, default)

        ## check validity of alternatives
        if (!any(is.na(alternatives)))
            for (i in alternatives)
                .check_value("alternatives", field, i)

        ## add field to meta data
        META <<- c(META, list(field))
        names(META)[length(META)] <<- name

        ## add (missing) fields to data entries
        DATA <<- lapply(DATA, `[[<-`, name, default)
    }

    .has_entry <-
    function(key)
        length(.get_entry_indices(key)) > 0

    .n_of_entries <-
    function()
        length(DATA)

    ## entry accessors
    .set_entry <-
    function(...)
    {
        ## check permissions
        if (!PERMISSIONS["set_entries"])
            stop("Setting of entries not allowed.", call. = FALSE)

        ## parameter handling
        l <- list(...)
        n <- names(l)
        index_fields <- .get_index_fields()
        fields <- .get_field_names()

        if (is.null(n)) {
            if (length(l) == length(fields)) {
                n <- fields
                names(l) <- n
            } else
                stop("Need either named arguments, or complete entry.")
        }


        .check_for_unknown_fields(n)

        ## check if there is at least one index field
        if (length(index_fields) < 1)
            stop("Need at least one index field before adding entries.",
                 call. = FALSE)

        ## check for mandatory fields
        mandatory_fields <- .get_mandatory_fields()
        missing_mandatory_fields <- !mandatory_fields %in% n
        if (any(missing_mandatory_fields))
            stop(paste("The following fields are mandatory, but missing:",
                       paste(mandatory_fields[missing_mandatory_fields],
                             collapse = ", ")), call. = FALSE)

        ## check for double entries
        if (length(.get_entry_indices(l[index_fields],
                                      stop_missing = FALSE,
                                      FUN = match_exact)) > 0)
                stop(paste("Entry already in registry."), call. = FALSE)

        ## check defaults and set values, if needed
        field_defaults    <- .get_field_defaults()
        default_fields    <- names(field_defaults)
        missing_fields    <- setdiff(default_fields, n)
        l[missing_fields] <- field_defaults[missing_fields]

        ## check field types, and apply field check function, if any.
        for (f in n)
            .check_value(f, .get_field(f), l[[f]])


        ## apply entry check function
        if (!is.null(validity_FUN))
            do.call(validity_FUN, list(l))

        ## add entry
        entry <- list(.make_entry(l))
        names(entry) <-
            paste(lapply(entry[[1]][seq_along(index_fields)], `[[`, 1),
                  collapse = "_")
        DATA <<- c(DATA, entry)
    }

    .get_entries <-
    function(...)
    {
        key <- list(...)
        if (length(key) > 0) {
            ind <- .get_entry_indices(key)
            if (length(ind) > 0)
                DATA[ind]
            else
                NULL
        } else DATA
    }

    .grep_entries <-
    function(pattern, ...)
    {
        pattern_in_entry <-
            function(x) any(sapply(x, function(i) is.character(i)
                                   && length(grep(pattern, i, ...) > 0)))
        ind <- sapply(DATA, pattern_in_entry)
        if (any(ind))
            DATA[ind]
        else
            NULL
    }

    .get_first_entry <-
    function(...)
        .get_entries(...)[[1]]

    .get_entry_names <-
    function()
    {
        if (length(DATA) < 1)
            character(0)
        else
            names(DATA)
    }

    .delete_entry <-
    function(...)
    {
        key <- list(...)

        ## check permissions
        if (!PERMISSIONS["delete_entries"])
            stop("Deletion of entries not allowed.", call. = FALSE)

        ## fetch entry index (this also checks if the entry exists)
        entry_index <- .get_entry_indices(key)

        ## check if it is unique
        if (length(entry_index) != 1)
            stop("Key specification must be unique.", call. = FALSE)

        ## check sealed entries
        if (entry_index %in% SEALED_ENTRIES)
            stop(paste("Deletion of entry not allowed."), call. = FALSE)

        ## delete it
        DATA[entry_index] <<- NULL
    }


    .modify_entry <-
    function(...)
    {
        ## check permissions
        if (!PERMISSIONS["modify_entries"])
            stop("Modifying of entries not allowed.", call. = FALSE)

        ## parameter handling
        l <- list(...)
        n <- names(l)

        .check_for_unknown_fields(n)

        ## determine entry key
        key  <- l[.get_index_fields()]

        ## fetch entry index (this also checks if the entry exists)
        entry_index <- .get_entry_indices(key)

        ## check if it is unique
        if (length(entry_index) != 1)
            stop("Key specification must be unique!", call. = FALSE)

        ## fetch entry
        entry <- DATA[[entry_index]]

        for(field in setdiff(n, .get_index_fields())) {
            ## check if field is modifiable
            field_entry <- .get_field(field)
            if (!field_entry$is_modifiable)
                stop(paste("Field", dQuote(field), "is not modifiable."),
                     call. = FALSE)

            ## check if entry and field are sealed
            if ((entry_index %in% SEALED_ENTRIES) &&
                (field %in% SEALED_FIELDS))
                stop(paste("Modification of field", dQuote(field),
                           "in this entry not allowed."), call. = FALSE)

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
            stop(paste("Field", dQuote(field), "not in registry."),
                 call. = FALSE)

        ret <- lapply(DATA, `[[`, field)
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
        PERMISSIONS["set_entries"] <<-
            PERMISSIONS["set_entries"] && set_entries
        PERMISSIONS["modify_entries"] <<-
            PERMISSIONS["modify_entries"] && modify_entries
        PERMISSIONS["delete_entries"] <<-
            PERMISSIONS["delete_entries"] && delete_entries
        PERMISSIONS["set_fields"] <<-
            PERMISSIONS["set_fields"] && set_fields
    }

    .seal_entries <-
    function()
    {
        SEALED_ENTRIES <<- seq_along(DATA)
        SEALED_FIELDS <<- .get_field_names()
    }

    .get_sealed_field_names <-
    function()
        SEALED_FIELDS

### CONSTRUCTOR
    ## return class
    structure(list(get_field = .FUNCall(.get_field),
                   get_fields = .FUNCall(.get_fields),
                   get_field_names = .FUNCall(.get_field_names),
                   set_field = .FUNCall(.set_field),

                   has_entry = .FUNCall(.has_entry),
                   get_entry = .FUNCall(.get_first_entry),
                   get_entries = .FUNCall(.get_entries),
                   get_entry_names = .FUNCall(.get_entry_names),
                   grep_entries = .FUNCall(.grep_entries),
                   set_entry = .FUNCall(.set_entry),
                   modify_entry = .FUNCall(.modify_entry),
                   delete_entry = .FUNCall(.delete_entry),
                   n_of_entries = .FUNCall(.n_of_entries),

                   get_field_entries = .FUNCall(.get_field_entries),

                   get_permissions = .FUNCall(.get_permissions),
                   restrict_permissions = .FUNCall(.restrict_permissions),
                   seal_entries = .FUNCall(.seal_entries),
                   get_sealed_field_names = .FUNCall(.get_sealed_field_names)
                   ),
              class = c(registry_class, "registry"))
}

