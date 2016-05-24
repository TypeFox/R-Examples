mock.conflicts <- structure(
    list(
        `http://local.crunch.io:8080/api/datasets/11dafed48fec42189511396b5f3c0c90/variables/b7112812955545b496570eb3745a5e6f/` = structure(
            list(
                conflicts = list(
                    structure(
                        list(
                            message = "Only in new dataset",
                            resolution = "Variable will be added with existing rows marked missing.",
                            resolved_by = NULL
                        ),
                        .Names = c("message", "resolution", "resolved_by")
                    )
                ),
                metadata = structure(
                    list(
                        references = structure(
                            list(
                                description = "", name = "mr_3",
                                format = structure(
                                    list(summary = structure(list(digits = 2), .Names = "digits")),
                                    .Names = "summary"
                                ),
                                dichotomous = TRUE, discarded = FALSE, alias = "mr_3", is_subvar = TRUE, header_order = NULL,
                                view = structure(
                                    list(show_counts = FALSE, column_width = NULL, include_missing = FALSE, show_numeric_values = FALSE),
                                    .Names = c("show_counts", "column_width", "include_missing", "show_numeric_values")
                                )
                            ),
                            .Names = c("description", "name", "format", "dichotomous", "discarded", "alias", "is_subvar", "header_order", "view")
                        ),
                        type = structure(
                            list(
                                class = "categorical",
                                categories = list(
                                    structure(
                                        list(numeric_value = 0, selected = FALSE, id = 1, name = "0.0", missing = FALSE),
                                        .Names = c("numeric_value", "selected", "id", "name", "missing")
                                    ),
                                    structure(
                                        list(numeric_value = 1, selected = TRUE, id = 2, name = "1.0", missing = FALSE),
                                        .Names = c("numeric_value", "selected", "id", "name", "missing")
                                    ),
                                    structure(
                                        list(numeric_value = NULL, selected = FALSE, id = -1, name = "No Data", missing = TRUE),
                                        .Names = c("numeric_value", "selected", "id", "name", "missing")
                                    )
                                ),
                                elements = list()
                            ),
                            .Names = c("class", "categories", "elements")
                        ),
                        derived = FALSE
                    ),
                    .Names = c("references", "type", "derived")
                )
            ),
            .Names = c("conflicts", "metadata")
        ),
        `http://local.crunch.io:8080/api/datasets/11dafed48fec42189511396b5f3c0c90/variables/bb79b8365526446bab7f6364ea111547/` = structure(
            list(
                conflicts = list(
                    structure(
                        list(
                            message = "Only in existing dataset",
                            resolution = "Additional rows will be marked missing.",
                            resolved_by = NULL
                        ),
                        .Names = c("message", "resolution", "resolved_by")
                    )
                ),
                metadata = structure(
                    list(
                        references = structure(
                            list(
                                description = "",
                                name = "mr_1",
                                format = structure(
                                    list(summary = structure(list(digits = 2), .Names = "digits")),
                                    .Names = "summary"
                                ),
                                dichotomous = TRUE,
                                discarded = FALSE,
                                alias = "mr_1",
                                is_subvar = TRUE,
                                header_order = NULL,
                                view = structure(
                                    list(show_counts = FALSE, column_width = NULL, include_missing = FALSE, show_numeric_values = FALSE),
                                    .Names = c("show_counts", "column_width", "include_missing", "show_numeric_values")
                                )
                            ),
                            .Names = c("description", "name", "format", "dichotomous", "discarded", "alias", "is_subvar", "header_order", "view")
                        ),
                        type = structure(
                            list(
                                class = "categorical",
                                categories = list(
                                    structure(
                                        list(numeric_value = 0, selected = FALSE, id = 1, name = "0.0", missing = FALSE),
                                        .Names = c("numeric_value", "selected", "id", "name", "missing")
                                    ),
                                    structure(
                                        list(numeric_value = 1, selected = TRUE, id = 2, name = "1.0", missing = FALSE),
                                        .Names = c("numeric_value", "selected", "id", "name", "missing")
                                    ),
                                    structure(
                                        list(numeric_value = NULL, selected = FALSE, id = -1, name = "No Data", missing = TRUE),
                                        .Names = c("numeric_value", "selected", "id", "name", "missing")
                                    )
                                ),
                                elements = list()
                            ),
                            .Names = c("class", "categories", "elements")
                        ),
                        derived = FALSE
                    ),
                    .Names = c("references", "type", "derived")
                )
            ),
            .Names = c("conflicts", "metadata")
        ),
        `http://local.crunch.io:8080/api/datasets/11dafed48fec42189511396b5f3c0c90/variables/0901dc49e8ab426bba64377148efcae2/` = structure(
            list(
                conflicts = list(
                    structure(
                        list(
                            message = "Subvariables didn't match",
                            resolution = "Union of subvariables will be used",
                            resolved_by = NULL
                        ),
                        .Names = c("message", "resolution", "resolved_by")
                    )
                ),
                metadata = structure(
                    list(
                        references = structure(
                            list(
                                description = "",
                                name = "MR",
                                format = structure(
                                    list(summary = structure(list(digits = 2), .Names = "digits")),
                                    .Names = "summary"
                                ),
                                dichotomous = TRUE,
                                discarded = FALSE,
                                alias = "MR",
                                is_subvar = NULL,
                                header_order = NULL,
                                view = structure(
                                    list(show_counts = FALSE, column_width = NULL, include_noneoftheabove = FALSE, include_missing = FALSE, show_numeric_values = FALSE),
                                    .Names = c("show_counts", "column_width", "include_noneoftheabove", "include_missing", "show_numeric_values")
                                )
                            ),
                            .Names = c("description", "name", "format", "dichotomous", "discarded", "alias", "is_subvar", "header_order", "view")
                        ),
                        type = structure(
                            list(
                                class = "categorical",
                                categories = list(
                                    structure(
                                        list(numeric_value = 0, selected = FALSE, id = 1, missing = FALSE, name = "0.0"),
                                        .Names = c("numeric_value", "selected", "id", "missing", "name")
                                    ),
                                    structure(
                                        list(numeric_value = 1, selected = TRUE, id = 2, missing = FALSE, name = "1.0"),
                                        .Names = c("numeric_value", "selected", "id", "missing", "name")
                                    ),
                                    structure(
                                        list(numeric_value = NULL, selected = FALSE, id = -1, missing = TRUE, name = "No Data"),
                                        .Names = c("numeric_value", "selected", "id", "missing", "name")
                                    )
                                ),
                                elements = list(),
                                matrix = c("a2b2e6f15a844ae59aecb753195ea528", "b7112812955545b496570eb3745a5e6f")
                            ),
                            .Names = c("class", "categories", "elements", "matrix")
                        ),
                        derived = TRUE
                    ),
                    .Names = c("references", "type", "derived")
                )
            ),
            .Names = c("conflicts", "metadata")
        )
    ),
    .Names = c(
        "http://local.crunch.io:8080/api/datasets/11dafed48fec42189511396b5f3c0c90/variables/b7112812955545b496570eb3745a5e6f/",
        "http://local.crunch.io:8080/api/datasets/11dafed48fec42189511396b5f3c0c90/variables/bb79b8365526446bab7f6364ea111547/",
        "http://local.crunch.io:8080/api/datasets/11dafed48fec42189511396b5f3c0c90/variables/0901dc49e8ab426bba64377148efcae2/"
    )
)
