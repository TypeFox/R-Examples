data("war1800")
f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | balanc + revis1 | balanc
m1 <- egame12(f1, data = war1800)

latexTable(m1)
latexTable(m1, digits = 8)
latexTable(m1, blankfill = "--")  ## Dashes in blank cells

\dontrun{
    latexTable(m1, file = "my_table.tex")  ## Write to file
}
