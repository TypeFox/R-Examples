sd_section(
    "Classes",
    "",
    c(
        "ColumnLinkedMatrix-class",
        "RowLinkedMatrix-class",
        "LinkedMatrix-class"
    )
)

sd_section(
    "Methods",
    "Methods for LinkedMatrix will be inherited by ColumnLinkedMatrix and RowLinkedMatrix.",
    c(
        "initialize,ColumnLinkedMatrix-method",
        "[<-,ColumnLinkedMatrix-method",
        "[,ColumnLinkedMatrix-method",
        "initialize,RowLinkedMatrix-method",
        "[<-,RowLinkedMatrix-method",
        "[,RowLinkedMatrix-method",
        "show,LinkedMatrix-method",
        "as.matrix.LinkedMatrix"
    )
)

sd_section(
    "Functions",
    "",
    c(
        "LinkedMatrix",
        "nNodes",
        "nodes",
        "index"
    )
)
