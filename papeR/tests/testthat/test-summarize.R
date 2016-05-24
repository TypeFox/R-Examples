library("papeR")
context("summarize functions")

if (require("nlme")) {
    ## Use dataset Orthodont
    data(Orthodont, package = "nlme")
    Ortho_small <<- Orthodont[Orthodont$Subject %in% c("M01", "M02", "F01", "F02"), ]

############################################################
## test old functions latex.table.fac / latex.table.cont
############################################################

test_that("latex.table.cont works", {
    expect_output(latex.table.cont(Orthodont),
                  paste0("tabular.*",
                         "& N &   & Mean & SD &   & Min & Q1 & Median & Q3 & Max",
                         ".*",
                         "distance & 108 &  & 24.02 & 2.93 &  & 16.50 & 22.00 & 23.75 & 26.00 & 31.50",
                         ".*",
                         "age & 108 &  & 11.00 & 2.25 &  & 8.00 & 9.00 & 11.00 & 13.00 & 14.00"))
    ## check that longtable isn't printed here
    expect_output(latex.table.cont(Orthodont), "(longtable){0}")
    ## but here
    expect_output(latex.table.cont(Orthodont, table = "longtable"), "longtable")
})

test_that("latex.table.fac works", {
    expect_output(latex.table.fac(Orthodont),
                  paste0("tabular.*",
                         "& Level &   & N & \\\\%",
                         ".*",
                         "Subject & M16 &  &   4 & 3.7",
                         ".*",
                         "Sex & Male &  &  64 & 59.3",
                         ".*",
                         "& Female &  &  44 & 40.7"))
    ## check that longtable isn't printed here
    expect_output(latex.table.fac(Orthodont), "(longtable){0}")
    ## but here
    expect_output(latex.table.fac(Orthodont, table = "longtable"), "longtable")
})

############################################################
## test variable labels in summaries
############################################################

test_that("variable labels work", {
    factor <- sapply(Orthodont, is.factor)
    for (type in c("numeric", "factor")) {
        data <- Ortho_small
        labels(data) <- c("Distance (mm)", "Age (yr)", "ID", "Sex")
        which <- if (type  == "numeric") { !factor } else { factor }

        ## summary with data set labels
        summary <- summarize(data, type = type, variable.labels = TRUE)
        expect_equivalent(summary[summary[, 1] != "", 1],
                          labels(data)[which],
                          info = type)

        ## summary with user specified labels
        usr_labels <- summarize(data, type = type,
                                variable.labels = c("a", "b", "c", "d"))
        expect_equivalent(usr_labels[usr_labels[, 1] != "", 1],
                          c("a", "b", "c", "d")[which],
                          info = type)

        ## grouped summary with data set labels
        which[4] <- FALSE
        grouped <- summarize(data, type = type, group = "Sex",
                             variable.labels = TRUE)
        expect_equivalent(grouped[grouped[, 1] != "", 1],
                          labels(data)[which],
                          info = type)
    }
})

test_that("grouped summaries work", {
    ## grouped summaries for numerics
    numeric <- summarize(Orthodont, type = "numeric", group = "Sex")
    expect_equivalent(numeric[, 2], rep(levels(Orthodont$Sex), 2))
    expect_equivalent(numeric$p.value[c(1,3)], c("<0.001", "1.000"))
    ## grouped summaries for factors
    factor <- summarize(Ortho_small, type = "factor", group = "Sex")
    expect_equivalent(factor[, 2], levels(Ortho_small$Subject))
    expect_equivalent(ncol(factor), 10)
    expect_equivalent(factor$p.value[1], "< 0.001")
})

test_that("print.summary works", {
    expect_output(print(summarize(Orthodont, type = "numeric")),
                  paste0("     N    Mean   SD    Min Q1 Median Q3  Max\n",
                         "1 distance 108   24.02 2.93   16.5 22  23.75 26 31.5\n",
                         "2      age 108   11.00 2.25    8.0  9  11.00 13 14.0"))
    expect_output(print(summarize(Orthodont, group = "Sex", type = "numeric")),
                  paste0("              Sex    N    Mean   SD    Min Q1 Median    Q3  Max   p.value\n",
                         "1 distance   Male   64   24.97 2.90   17.0 23  24.75 26.50 31.5    <0.001\n",
                         "2          Female   44   22.65 2.40   16.5 21  22.75 24.25 28.0          \n",
                         "3      age   Male   64   11.00 2.25    8.0  9  11.00 13.00 14.0     1.000\n",
                         "4          Female   44   11.00 2.26    8.0  9  11.00 13.00 14.0          "))
    expect_output(print(summarize(Orthodont, type = "factor")),
                  paste0("            Level    N    %\n",
                         "1  Subject    M16    4  3.7\n",
                         "2             M05    4  3.7\n",
                         "3             M02    4  3.7\n",
                         ".*",
                         "28     Sex   Male   64 59.3\n",
                         "29         Female   44 40.7"))
    expect_output(print(summarize(Ortho_small, group = "Sex", type = "factor")),
                  paste0("                  Sex: Male        Sex: Female               \n",
                         "          Level           N    %             N    %   p.value\n",
                         "1 Subject   M02           4 50.0             0  0.0   < 0.001\n",
                         "2           M01           4 50.0             0  0.0          \n",
                         "3           F01           0  0.0             4 50.0          \n",
                         "4           F02           0  0.0             4 50.0          "))

})

test_that("caption works", {
    ## via call arguments
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test"),
                        floating = TRUE),
                  ".*\\caption\\{Test\\}.*")
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test",
                               label= "tab:Test"),
                        floating = TRUE),
                  ".*\\caption\\{Test\\}.*\\label\\{tab:Test\\}")

    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test"),
                        tabular.environment = "longtable", floating = TRUE),
                  ".*\\caption\\{Test\\}.*")
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test",
                               label= "tab:Test"),
                        tabular.environment = "longtable", floating = TRUE),
                  ".*\\caption\\{Test\\}.*\\label\\{tab:Test\\}")

    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test"),
                        tabular.environment = "longtable"),
                  ".*\\caption\\{Test\\}.*")
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test",
                               label= "tab:Test"),
                        tabular.environment = "longtable"),
                  ".*\\caption\\{Test\\}.*\\label\\{tab:Test\\}")

    ## requires capt-of
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test")),
                  paste0(".*begin\\{minipage\\}\\{.*linewidth\\}\n",
                         ".*captionof\\{table\\}\\{Test\\}\n",
                         ".*",
                         "end\\{minipage\\}"))
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"), caption= "Test",
                               label= "tab:Test")),
                  paste0(".*begin\\{minipage\\}\\{.*linewidth\\}\n",
                         ".*captionof\\{table\\}\\{Test\\}\n",
                         ".*label\\{tab:Test\\}\n",
                         ".*",
                         "end\\{minipage\\}"))
    ## additionally test if this also works via options
})

test_that("endhead is included if necessary", {
    ## via call arguments
    expect_output(print(xtable(summarize(Orthodont, type = "numeric")),
                        tabular.environment = "longtable"),
                  ".*cmidrule\\{7-11\\}\n.*endhead\ndistance.*")
    ## via options
    options(xtable.tabular.environment = "longtable")
    expect_output(print(xtable(summarize(Orthodont, type = "numeric"))),
                  ".*cmidrule\\{7-11\\}\n.*endhead\ndistance.*")
    options(xtable.tabular.environment = NULL)
})

test_that("xtable works for summarize_factor with groups", {
    grouped <- summarize(Ortho_small, type = "factor", group = "Sex")
    expect_output(print(xtable(grouped)),
                  paste(".* & & &\\\\multicolumn\\{2\\}\\{c\\}\\{Sex: Male\\} &",
                        "& \\\\multicolumn\\{2\\}\\{c\\}\\{Sex: Female\\} &  &.*"))
    grouped <- summarize(Ortho_small, type = "factor", group = "Sex", test = FALSE)
    expect_output(print(xtable(grouped)),
                  paste(".* & & &\\\\multicolumn\\{2\\}\\{c\\}\\{Sex: Male\\} &",
                        "& \\\\multicolumn\\{2\\}\\{c\\}\\{Sex: Female\\}.*"))
})

test_that("include.rownames is ignored", {
    tab <- summarize(Orthodont, type = "numeric")
    expect_output(print(xtable(tab), include.rownames = FALSE),
                  ".*\n distance .*\n  age .*")
    expect_warning(print(xtable(tab), include.rownames = TRUE))
    expect_output(print(xtable(tab), include.rownames = TRUE),
                  ".*\n distance .*\n  age .*")
    expect_output(print(xtable(tab)),
                  ".*\n distance .*\n  age .*")
})

}


