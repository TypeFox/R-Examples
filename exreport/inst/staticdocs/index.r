sd_section("Experimental Data Manipulation",
           "This functions manipulate the experiment objects, implementing
           common data manipulation operations, such as subseting and joining.
           In addition, there are a number of available functions to sanitize
           and preprecessing the experimental data such as looking for duplicated
           entries, reducing parameters and instantiating the methods with the
           parameters.",
           c(
             "expCreate",
             "expCreateFromTable",
             
             "expCombine",
             "expConcat",
             "expExtend",
             "expExtract",
             "expGetDuplicated",
             "expInstantiate",
             "expReduce",
             "expRemoveDuplicated",
             "expRename",
             "expReorder",
             "expSubset"
           )
)

sd_section("Statistical Tests",
           "This functions implement several statistical test to compare the
           methods of the experiment.",
           c(
             "testPaired",
             "testMultipleControl",
             "testMultiplePairwise"
           )
)

sd_section("Tabular Data Generation",
           "This functions generate tables summarizing the information of an
           experiment or a text to be printed to pdf or web reports.",
           c(
             "tabularExpSummary",
             "tabularTestPairwise",
             "tabularTestSummary"
           )
)

sd_section("Graphical Plots",
           "This functions generate plots summarizing the information of an
           experiment or a text to be printed to pdf or web reports.",
           c(
             "plotCumulativeRank",
             "plotExpSummary",
             "plotRankDistribution"
           )
)

sd_section("Rendering Reports",
           "These are the main functions used to generate and render the reports",
           c(
             "exreport",
             "exreportAdd",
             "exreportRender"
           )
)

sd_section("Problems",
           "Example problems for the examples and documentations",
           c(
             "wekaExperiment"
           )
)