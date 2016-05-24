# Create file that declares "global" variables in figure*() code
# (actually just variables that codetools cannot find)
if(getRversion() >= "2.15.1") {
    utils::globalVariables(unique(c(
        "supp", "grayify", "year", "displ", "hwy", "cyl", "YHOO", "GlaucomaM",
        "myNode", "soil", "NHANES", "Category", "Count", "Var1", "Var2",
        "Freq", "grd", "DEPTH", "DEPTH", "DEPTH", "DEPTH", "gear", "gear",
        "lowTideDate", "lowTideHour", "mainHours", "phases", "gear", "gear",
        "temperature", "trans", "cyl", "disp", "trans", "disp", "disp",
        "disp", "limit", "category", "temperature", "disp", "gear", "disp",
        "trans", "trans", "disp", "disp", "trans", "gear", "disp", "trans",
        "disp", "disp")))
    # Data sets from 'datasets' that cannot be import()ed because
    # 'datasets' does not export anything
    utils::globalVariables(c("datasets", "LifeCycleSavings", "Nile",
                             "OrchardSprays", "Titanic", "ToothGrowth",
                             "USArrests", "USJudgeRatings", "VADeaths",
                             "faithful", "iris", "mtcars", "nhtemp",
                             "pressure", "trees", "volcano"))
}

