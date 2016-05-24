run.integration.tests <- Sys.getenv("INTEGRATION") == "TRUE"
Sys.setlocale("LC_COLLATE", "C") ## What CRAN does

skip_on_jenkins <- function (...) {
    if (nchar(Sys.getenv("JENKINS_HOME"))) {
        skip(...)
    }
}

skip_locally <- function (...) {
    if (nchar(Sys.getenv("JENKINS_HOME")) == 0) {
        skip(...)
    }
}

set.seed(666)

# httpcache::startLog("") ## prints to stdout

fromJSON <- jsonlite::fromJSON
loadLogfile <- httpcache::loadLogfile
cacheLogSummary <- httpcache::cacheLogSummary
requestLogSummary <- httpcache::requestLogSummary
uncached <- httpcache::uncached

envOrOption <- function (opt) {
    ## .Rprofile options are like "test.api", while env vars are "R_TEST_API"
    envvar.name <- paste0("R_", toupper(gsub(".", "_", opt, fixed=TRUE)))
    envvar <- Sys.getenv(envvar.name)
    if (nchar(envvar)) {
        ## Let environment variable override .Rprofile, if defined
        return(envvar)
    } else {
        return(getOption(opt))
    }
}

## .onAttach stuff, for testthat to work right
options(
    crunch.api=envOrOption("test.api"),
    warn=1,
    crunch.debug=FALSE,
    digits.secs=3,
    crunch.timeout=15,
    httpcache.on=TRUE,
    # httpcache.log="",
    crunch.namekey.dataset="alias",
    crunch.namekey.array="name",
    crunch.email=envOrOption("test.user"),
    crunch.pw=envOrOption("test.pw")
)
set_config(crunchConfig())

## Test serialize and deserialize
cereal <- function (x) fromJSON(toJSON(x), simplifyVector=FALSE)

newDatasetFromFixture <- function (filename) {
    ## Grab csv and json from "dataset-fixtures" and make a dataset
    m <- fromJSON(file.path("dataset-fixtures", paste0(filename, ".json")),
        simplifyVector=FALSE)
    return(suppressMessages(createWithMetadataAndFile(m,
        file.path("dataset-fixtures", paste0(filename, ".csv")))))
}

releaseAndReload <- function (dataset) {
    .releaseDataset(dataset)
    return(refresh(dataset))
}

## Data frames to make datasets with
df <- data.frame(v1=c(rep(NA_real_, 5), rnorm(15)),
                 v2=c(letters[1:15], rep(NA_character_, 5)),
                 v3=8:27,
                 v4=as.factor(LETTERS[2:3]),
                 v5=as.Date(0:19, origin="1955-11-05"),
                 v6=TRUE,
                 stringsAsFactors=FALSE)

mrdf <- data.frame(mr_1=c(1,0,1,NA_real_),
                   mr_2=c(0,0,1,NA_real_),
                   mr_3=c(0,0,1,NA_real_),
                   v4=as.factor(LETTERS[2:3]),
                   stringsAsFactors=FALSE)

testfile.csv <- "fake.csv"
testfile.df <- read.csv(testfile.csv)

mrdf.setup <- function (dataset, pattern="mr_", name=ifelse(is.null(selections),
                        "CA", "MR"), selections=NULL) {
    cast.these <- grep(pattern, names(dataset))
    dataset[cast.these] <- lapply(dataset[cast.these],
        castVariable, "categorical")
    if (is.null(selections)) {
        dataset[[name]] <- makeArray(pattern=pattern, dataset=dataset,
            name=name)
    } else {
        dataset[[name]] <- makeMR(pattern=pattern, dataset=dataset, name=name,
            selections=selections)
    }
    return(dataset)
}

validImport <- function (ds) {
    ## Pull out common tests that "df" was imported correctly
    expect_true(is.dataset(ds))
    expect_identical(description(ds), "")
    expect_identical(names(df), names(ds))
    expect_identical(dim(ds), dim(df))
    expect_true(is.Numeric(ds[["v1"]]))
    expect_true(is.Text(ds[["v2"]]))
    expect_identical(name(ds$v2), "v2")
    expect_true(is.Numeric(ds[["v3"]]))
    expect_identical(description(ds$v3), "")
    expect_equivalent(as.array(crtabs(mean(v3) ~ v4, data=ds)),
        tapply(df$v3, df$v4, mean, na.rm=TRUE))
    expect_equivalent(as.vector(ds$v3), df$v3)
    expect_true(is.Categorical(ds[["v4"]]))
    expect_equivalent(as.array(crtabs(~ v4, data=ds)),
        array(c(10, 10), dim=2L, dimnames=list(v4=c("B", "C"))))
    expect_true(all(levels(df$v4) %in% names(categories(ds$v4))))
    expect_identical(categories(ds$v4), categories(refresh(ds$v4)))
    expect_identical(ds$v4, refresh(ds$v4))
    expect_equivalent(as.vector(ds$v4), df$v4)
    expect_true(is.Datetime(ds$v5))
    expect_true(is.Categorical(ds$v6))
    expect_identical(showVariableOrder(ordering(ds)), names(variables(ds)))
}

validApidocsImport <- function (ds) {
    expect_true(is.dataset(ds))
    expect_identical(dim(ds), c(20L, 9L))
    expect_identical(names(ds),
        c("allpets", "q1", "petloc", "ndogs", "ndogs_a", "ndogs_b", "q3",
        "country", "wave"))

}

## Global teardown
bye <- new.env()
if (run.integration.tests) {
    with(test.authentication, {
        datasets.start <- urls(datasetCatalog())
        users.start <- urls(getUserCatalog())
        projects.start <- urls(session()$projects)
    })
}
reg.finalizer(bye,
    function (x) {
        if (run.integration.tests) {
            with(test.authentication, {
                datasets.end <- urls(datasetCatalog())
                leftovers <- setdiff(datasets.end, datasets.start)
                if (length(leftovers)) {
                    stop(length(leftovers),
                        " dataset(s) created and not destroyed: ",
                        serialPaste(dQuote(names(datasetCatalog()[leftovers]))),
                        call.=FALSE)
                }
                users.end <- urls(getUserCatalog())
                leftovers <- setdiff(users.end, users.start)
                if (length(leftovers)) {
                    stop(length(leftovers),
                        " users(s) created and not destroyed: ",
                        serialPaste(dQuote(names(getUserCatalog()[leftovers]))),
                        call.=FALSE)
                }
                projects.end <- urls(session()$projects)
                leftovers <- setdiff(projects.end, projects.start)
                if (length(leftovers)) {
                    stop(length(leftovers),
                        " projects(s) created and not destroyed: ",
                        serialPaste(dQuote(names(session()$projects[leftovers]))),
                        call.=FALSE)
                }
            })
        }
    },
    onexit=TRUE)
