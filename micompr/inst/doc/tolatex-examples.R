## ----micdata, results = 'asis', warning = FALSE--------------------------
library(micompr)
mic <- micomp(4, 0.95,
             list(list(name = "NLvsJOK", grpout = pphpc_ok),
                  list(name = "NLvsJNS", grpout = pphpc_noshuff),
                  list(name = "NLvsJDIF", grpout = pphpc_diff)),
             concat = TRUE)

## ----table01, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        caption = "Default table.")

## ----table02, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        caption = "Booktabs.", booktabs = T)

## ----table03, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        booktabs = T, labels_cmp_show = F,
        caption = "No comparison label.")

## ----table04, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        booktabs = T, labels_col_show = F,
        caption = "No data label.")

## ----table05, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        booktabs = T, labels_cmp_show = F, labels_col_show = F,
        caption = "No data and comparison labels.")

## ----table06, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        tag_comp = "Comparisons", tag_data = "What?",
        tag_outputs = "Outs.",
        data_labels = c("No. PCs", "MANOVA", "$t$-test",
                        "Mann-Whitney", "PC1 vs PC2"),
        caption = "Alternative header tags and data labels.")

## ----table07, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        label_row_show = F,
        caption = "Do not show outputs tag.")

## ----table08, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        booktabs = T, labels_cmp_show = F,
        labels_col_show = F, label_row_show = F,
        caption = paste0("No data and comparison labels and ",
                         "no outputs tag, with booktabs."))

## ----table09, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        data_show = c("parp-1", "parp-2", "sep",
                      "aparp-1", "aparp-2", "sep",
                      "varexp-1", "varexp-2"),
        caption = "Different types of data, with separators.")

## ----table10, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        booktabs = T,
        data_show = c("parp-1", "parp-2", "sep",
                      "aparp-1", "aparp-2", "sep",
                      "varexp-1", "varexp-2"),
        data_labels = c("$t$-test 1", "$t$-test 2",
                      "$t$-test 1 (wb)", "$t$-test 2 (wb)",
                      "Var 1", "Var 2"),
        caption = paste0("Different types of data, booktabs, ",
                         "custom data labels."))

## ----table11, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        orientation = F,
        data_labels = c("NoPCs", "MNV", "$t$", "MW", NA),
        scoreplot_before =
          "\\raisebox{-.5\\height}{\\resizebox {0.7cm} {0.7cm} {",
        caption = paste0("Transposed table with score plots and ",
                         "NA in one of the data labels (such ",
                         "that a default should be used)."))

## ----table12, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        orientation = F,
        booktabs = T,
        data_show = c("npcs-1", "mnvp-1", "parp-1", "nparp-1"),
        data_labels = c("NoPCs", "MNV", "$t$", "MW"),
        caption = paste0("Transposed table, without score ",
                         "plots, with booktabs."))

## ----table13, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        orientation = F,
        booktabs = T,
        pvalf_params = list(minval = 1e-6, na_str = "$\\times$"),
        data_show = c("npcs-1", "mnvp-1", "parp-1", "nparp-1"),
        data_labels = c("NoPCs", "MNV", "$t$", "MW"),
        labels_cmp_show = F, labels_col_show = F,
        label_row_show = F,
        caption = paste0("Transposed table: without score ",
                         "plots, with booktabs, custom ",
                         "p-value parameters."))

## ----table14, results = 'asis', warning = FALSE--------------------------
toLatex(mic,
        orientation = F,
        booktabs = T,
        data_show = c("parp-1", "parp-2", "sep",
                      "aparp-1", "aparp-2", "sep",
                      "varexp-1", "varexp-2"),
        data_labels = c("$t_1$", "$t_2$",
                      "$t_1\\ast$", "$t_2\\ast$",
                      "$V_1$", "$V_2$"),
        caption = paste0("Transposed table, different types ",
                         "of data, booktabs, ",
                         "custom data labels."))

## ----table15, results = 'asis', warning = FALSE--------------------------
toLatex(mic[[1, 1]],
        orientation = F,
        labels_cmp_show = F,
        label_row_show = F,
        booktabs = T,
        data_show = c("npcs-1", "mnvp-1", "parp-1", "nparp-1", "scoreplot"),
        data_labels = c("NoPCs", "MNV", "$t$", "MW", "Scores"),
        caption = paste0("Table with a single cmpoutput object."))

