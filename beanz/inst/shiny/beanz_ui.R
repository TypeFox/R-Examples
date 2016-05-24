##
##  UI definitions for BEANZ project
##

DATA.FORMAT <- c("Summary treatment effect data", "Subject level raw data");
DISPLAY.BY  <- c('By results', 'By model');
ALL.MODELS  <- c("No subgroup effect",
                 "Full stratification",
                 "Simple regression",
                 "Basic shrinkage",
                 "Regression and shrinkage",
                 "Dixon and Simon",
                 "Extended Dixon and Simon"
                 );
STAN.NAME   <- c("nse", "fs", "sr", "bs", "srs", "ds", "eds");
RESP.TYPE   <- c("Continuous",
                 "Binary with treatment effect measured by log odds ratio",
                 "Time to event with treatment effect measured by log hazard ratio");
SUB.HEAD    <- c("Estimate", "Variance");


##-------------------------------------------------------------
##           UI DEFINITIONS
##-------------------------------------------------------------

##show different type of messages
msg.box <- function(contents, type="info") {
    switch(type,
           info    = cls <- "cinfo",
           warning = cls <- "cwarning",
           success = cls <- "csuccess",
           error   = cls <- "cerror");
    rst <- '<div class="';
    rst <- paste(rst, cls, '">');
    rst <- paste(rst, contents);
    rst <- paste(rst, "</div>");
    HTML(rst);
}

##chk number of subgroups to determine plot or not
gen.chk.nsub <- function(x, max.s, A) {
    if (x > max.s) {
        msg.box("There are too many subgroup comparisons.
                     Density plots are not generated.", "warning")
    } else if (x < 2) {
        msg.box("There is only one subgroup selected.
                     Density plots are not generated.", "warning")
    } else {
        A
    }
}

##generate analysis result panels
gen.rst.tabs <- reactive({
    mrst      <- ana.rst();
    n.sel     <- length(get.subgrp.selected());
    max.s     <- input$displaynsub;

    if (is.null(mrst))
        return(NULL);

    rst <- list(widths=c(3,9));
    for (i in 1:length(mrst)) {
        rst[[length(rst)+1]] <- names(mrst)[i];
        rst[[length(rst)+1]] <- tabPanel("Effects",
                                         tabsetPanel(position = "above", type="pills",
                                                     tabPanel("Effects Table",
                                                              DT::dataTableOutput(paste("tblrst",i,sep=""))),
                                                     tabPanel("Density",
                                                              plotOutput(paste("plotrst",i,sep=""))),
                                                     tabPanel("Forest plot",
                                                              plotOutput(paste("plotforest",i,sep=""))),
                                                     align="center"
                                                     )
                                         );
        rst[[length(rst)+1]] <- tabPanel("Comparisons",
                                         tabsetPanel(position = "above", type="pills",
                                                     tabPanel("Effects Comparison Table",
                                                              DT::dataTableOutput(paste("tblcomp",i,sep=""))),
                                                     tabPanel("Density",
                                                              lapply(n.sel,
                                                                     gen.chk.nsub,
                                                                     max.s,
                                                                     plotOutput(paste("plotcomp",i,sep="")))
                                                              ),
                                                     tabPanel("Forest plot",
                                                              lapply(n.sel,
                                                                     gen.chk.nsub,
                                                                     max.s,
                                                                     plotOutput(paste("plotforestcomp",i,sep=""))
                                                                     )),
                                                     align="center")
                                         );
        rst[[length(rst)+1]] <- tabPanel("STAN Diagnosis",
                                         tabsetPanel(position = "above", type="pills",
                                                     tabPanel("DIC",
                                                              htmlOutput(paste("txtdic",i,sep=""))),
                                                     tabPanel("Output",
                                                              verbatimTextOutput(paste("txtrst",i,sep=""))),
                                                     tabPanel("Trace Plot",
                                                              plotOutput(paste("plottrace",i,sep=""))
                                                              ))
                                         );
    }
    rst
})

##genereate individual tabpanel for each model
get.stab <- function(mrst, lab, f) {
    do.call(tabsetPanel,
            c(list(position = "above", type="pills"),
              lapply(1:length(mrst),
                     function(i) {
                         tabPanel(names(mrst)[i],
                                  f(paste(lab,i,sep="")))
                     })
              ));
}

##across by models
gen.rst.tabs.mdl <- reactive({
    mrst  <- ana.rst();
    n.sel <- length(get.subgrp.selected());
    max.s <- input$displaynsub;

    if (is.null(mrst))
        return(NULL);

    pan.table      <- get.stab(mrst, "tblrst", DT::dataTableOutput);
    pan.den        <- get.stab(mrst, "plotrst", plotOutput);
    pan.forest     <- get.stab(mrst, "plotforest", plotOutput);
    pan.tablecomp  <- get.stab(mrst, "tblcomp", DT::dataTableOutput);
    pan.stan       <- get.stab(mrst, "txtrst", verbatimTextOutput);
    pan.trace      <- get.stab(mrst, "plottrace", plotOutput);
    pan.dic        <- get.stab(mrst, "txtdic", htmlOutput);

    if (n.sel > max.s) {
        pan.dencomp <- pan.forestcomp <- msg.box("There are too many subgroup comparisons.
                                                   Plots are not generated.", "warning")
    } else if (n.sel < 2) {
        pan.dencomp <- pan.forestcomp <- msg.box("There is only one subgroup selected.
                                                  Plot is not generated.", "warning")

    } else {
        pan.dencomp    <- get.stab(mrst, "plotcomp", plotOutput);
        pan.forestcomp <- get.stab(mrst, "plotforestcomp", plotOutput);
    }


    rst <- list(widths=c(2,10),
                "Effect",
                tabPanel("Table", pan.table),
                tabPanel("Density", pan.den),
                tabPanel("Forest plot", pan.forest),
                "Comparison",
                tabPanel("Table", pan.tablecomp),
                tabPanel("Density", pan.dencomp),
                tabPanel("Forest plot", pan.forestcomp),
                "STAN Diagnosis",
                tabPanel("DIC", pan.dic),
                tabPanel("Output", pan.stan),
                tabPanel("Trace Plot", pan.trace)
                );
})


##subgroup panel
panel.subgroup <- function() {
    list(
        wellPanel(
            h4(textOutput("txtSel")),
            msg.box("Specify variables for defining subgroups.
                     If a covariate has more than two levels in the data,
                     please specify whether the covariate is nomial or ordinal.
                     Ordinal covariates are assumed to have a linear trend across
                     levels in its effect on the outcome."),
            fluidRow(
                conditionalPanel(paste("input.dataformat=='", DATA.FORMAT[2], "'",sep=""),
                                 column(3,
                                        h6("Treatment"),
                                        uiOutput("uiTrt"),
                                        h6("Response"),
                                        uiOutput("uiResp"),
                                        conditionalPanel(paste("input.resptype=='", RESP.TYPE[3], "'",sep=""),
                                                         h6("Censoring"),
                                                         uiOutput("uiCensor"))
                                        ),
                                 column(2,
                                        h6("Covariates"),
                                        uiOutput("uiCov")
                                        ),
                                 column(3,
                                        h6("Type of response"),
                                        radioButtons('resptype', '', RESP.TYPE)
                                        )
                                 ),
                conditionalPanel(paste("input.dataformat=='", DATA.FORMAT[1], "'",sep=""),
                                 column(3,
                                        h6("Treatment Effect Estimates"),
                                        uiOutput("uisubE"),
                                        h6("Treatment Effect Variance"),
                                        uiOutput("uisubVar")),
                                 column(2,
                                        h6("Covariates"),
                                        uiOutput("uisubCov"))),
                column(3, uiOutput("uiNominal"))
            ),
            fluidRow(
                column(9,
                       actionButton("btnSub", "Get Subgroups", styleclass = "success"),
                       align="left")
            ),

            htmlOutput("uiValid")
        ),
        ##subgroup
        uiOutput("uiSubgrp")
    )
}

panel.priors <- function() {

    ##no subgroup effect label
    n.lab1 <- sprintf("<span class='hline2'>
                           %s
                          </span>", ALL.MODELS[1]);


    f.each <- function(inx, ...) {

        n.chk <- paste("chkM", inx, sep="");
        n.par <- paste("parM", inx, sep="");
        n.lab <- sprintf("<span class='hline2'>
                           %s
                          </span>", ALL.MODELS[inx])

        withTags(div(class="content",
                     checkboxInput(inputId = n.chk,
                                   label = HTML(n.lab),
                                   value=TRUE),
                     div(class='intro', id=n.par, ...))
                 )
    }

    fluidRow(column(4,
                    withTags(div(class="content",
                                 HTML(n.lab1),
                                 div(class='intro', id="parM1",
                                     h6("Prior variance of the group mean"),
                                     sliderInput(inputId = "m1vtau", label = "", value = 1000, min = 100, max = 10000, step=100)
                                     ))
                             ),
                    f.each(2,
                          h6("Prior variance of the group mean"),
                           sliderInput(inputId = "m2vtau", label = "", value = 1000, min = 100, max = 10000, step=100)
                           ),
                    f.each(3,
                           h6("Prior variance of the intercept"),
                           sliderInput(inputId = "m3vtau", label = "", value = 1000, min = 100, max = 10000, step=100),
                           h6("Prior variance of regression coefficients"),
                           sliderInput(inputId = "m3vgamma", label = "", value = 1000, min = 100, max = 10000, step=100)
                           )
                    ),
             column(4,
                    f.each(4,
                           h6("Prior variance of the group mean"),
                           sliderInput(inputId = "m4vtau", label = "", value = 1000, min = 100, max = 10000, step=100),
                           h6("Prior variance of the shrinkage parameter"),
                           sliderInput(inputId = "m4vw", label = "", value = 100, min = 0.1, max = 500, step=0.5)
                           ),
                    f.each(5,
                           h6("Prior variance of the intercept"),
                           sliderInput(inputId = "m5vtau", label = "", value = 1000, min = 100, max = 10000, step=100),
                           h6("Prior variance of regression coefficients"),
                           sliderInput(inputId = "m5vgamma", label = "", value = 1000, min = 100, max = 10000, step=100),
                           h6("Prior variance of the shrinkage parameter"),
                           sliderInput(inputId = "m5vw", label = "", value = 100, min = 0.1, max = 500, step=0.5)
                          )
                    ),
             column(4,
                    f.each(6,
                           h6("Prior variance of the group mean"),
                           sliderInput(inputId = "m6vtau", label = "", value = 1000, min = 100, max = 10000, step=100),
                           h6("Prior variance of the shrinkage parameter"),
                           sliderInput(inputId = "m6vw", label = "", value = 100, min = 0.1, max = 500, step=0.5)
                           ),
                    f.each(7,
                           h6("Prior variance of the group mean"),
                           sliderInput(inputId = "m7vtau", label = "", value = 1000, min = 100, max = 10000, step=100),
                           h6("Prior variance of the shrinkage parameter"),
                           sliderInput(inputId = "m7vw", label = "", value = 100, min = 0.1, max = 500, step=0.5)
                           )
                    )
             );
}


##tabset for user log in
tab.login <- function(){
    div(class="well",
        textInput("inUname", "User Name", value="guest"),
        textInput("inPwd", "Password", value="guest"),
        br(),
        actionButton("btnLogin", "Login", styleclass="info"),
        style='width:200px; margin:0 auto; min-height:200px; border-width:0')
}


##entrance page
tab.entrance <- function() {
    HTML('<div class="entrance">
          <img alt="compact" src="compact.svg" data-img="svg"></img>
                    <h4>Light Version</h4>
                        <p> For non-statisticians,
                            options optimized </p>
                            <button id="btnLight" type="button" class="btn action-button btn-info shiny-bound-input">
                                Visit
                            </button>
        </div>
        <div class="entrance">
          <img alt="compact" src="fullversion.svg" data-img="svg"></img>
                    <h4>Full Version</h4>
                        <p> For statisticians, options available </p>
                            <button id="btnFull" type="button" class="btn action-button btn-info shiny-bound-input">
                                Visit
                            </button>
        </div>')

}

##tabset for start page
tab.start <- function() {
    tabPanel("Start",
             navlistPanel(widths=c(2,10),
                          tabPanel("What does BEANZ do?",
                                   wellPanel(includeHTML("www/beanz_do.html"))),
                          tabPanel("What does BEANZ need?",
                                   wellPanel(includeHTML("www/beanz_need.html"))),
                          tabPanel("What does BEANZ provide?",
                                   wellPanel(includeHTML("www/beanz_rst.html"))))
             );
}


##tabset for data uploading
tab.upload <- function(){
    tabPanel("Upload Data",
             fluidPage(
                 ##msg.box('Please upload data file on this page.
                 ##         Note at this stage, the software only considers discrete baseline covariates.
                 ##         Click to download a <a href="example_rawdata.txt" download>raw data</a>
                 ##         or <a href="example_subgroupinfo.txt" download>subgroup treatment effect </a> example file.'),
                 wellPanel(
                     h4("Upload data"),
                     msg.box("Please upload data file."),
                     fluidRow(
                         column(3,
                                h6("Choose File"),
                                fileInput(inputId = 'userdata', label = '',
                                          accept=c('text/csv','text/comma-separated-values,text/plain'))),
                         column(3,
                                h6("Data Format"),
                                radioButtons('dataformat', '', DATA.FORMAT, selected = DATA.FORMAT[2])
                                ),
                         column(2, h6("Separator"),
                                radioButtons('sep', '',
                                             c(Comma=',',Semicolon=';',Tab='\t',Space=' '),
                                             ' ')),
                         column(2, h6("Quote"),
                                radioButtons('quote', '',
                                             c(None='','Double Quote'='"','Single Quote'="'"),
                                             selected = '')),
                         column(2, h6("Other"),
                                checkboxInput(inputId='header', label='Header', value=TRUE),
                                checkboxInput(inputId="displaydata", label = "Show Data", value = TRUE))
                     )),
                 wellPanel(
                     h4("Try An Example"),
                     includeHTML("www/beanz_solvd.html"),
                     actionButton("btnExample", "Try it", styleclass = "success")
                 ),
                 wellPanel(h4("Review data"), DT::dataTableOutput("uiData"))
                 )
             )
}

##select models
tab.model <- function() {
    tabPanel("Statistical Models",
             fluidPage(
                 wellPanel(
                     msg.box("Select statistical models for subgroup analysis."),
                     fluidRow(
                         column(3,
                                h4("Models"),
                                checkboxGroupInput("selmdl","",
                                                   choices = ALL.MODELS,
                                                   selected = ALL.MODELS[1:length(ALL.MODELS)])),
                         column(9,
                                h4("Statistical Details"),
                                htmlOutput("uiMdlInst")))))
             );
}

##mcmc configuration
tab.mcmc <- function() {
    tabPanel("Configuration",
             wellPanel(
                 h4("Statistical Models and Priors"),
                 msg.box("<p>Select statistical models and specify priors for model parameters.</p>
                          <p>Details of the models can be found
                          in the <a href='paper_beanz.pdf'>software manual</a>.</p>"),
                 panel.priors()
             ),
             wellPanel(
                 h4("MCMC Paramters"),
                 msg.box("Specify parameters for Bayesian posterior sampling."),
                 fluidRow(
                     column(4,
                            h6("Number of iterations"),
                            sliderInput(inputId = "mcmciter", label = "",
                                        value = 4000, min = 100, max = 20000, step=100),
                            h6("Number of burn-in"),
                            sliderInput("mcmcburnin", label = "", value=2000, min=100, max=20000, step=100),
                            h6("Number of thinning"),
                            sliderInput(inputId = "mcmcthin", label = "", value=2, min=1, max=50, step=1)),
                     column(3,
                            h6("Random seed"),
                            numericInput(inputId="mcmcseed", label="", value=0, min=0),
                            h6("Algorithm"),
                            radioButtons('mcmcalg', '', c('NUTS', 'HMC', "Fixed_param"))
                            ),
                     column(3,
                            h6("Uncertainty of log variance"),
                            sliderInput(inputId = "rangelogvar", label = "",
                                        min = -1, max = 1, value = c(-0.000001, 0.0000001), step = 0.01)
                            )
                 )),
             wellPanel(
                 h4("Display Parameters"),
                 fluidRow(
                     column(4,
                            h6("Cut off value for treatment effects"),
                            numericInput(inputId="displaycut", label="", value=0),
                            h6("Cut off value for subgroup comparison"),
                            numericInput(inputId="displaycutcomp", label="", value=0),
                            h6("Digits"),
                            numericInput(inputId="displaydigit", label="", value=3, min=1, step=1)
                            ),
                     column(4,
                            h6("Maximum subgroups for comparison plots"),
                            numericInput(inputId="displaynsub", label="", value=5, step=1),
                            h6("Organize results"),
                            radioButtons('displayby', '', DISPLAY.BY)
                            ),
                     column(4,
                            h6("Transformation"),
                            checkboxInput(inputId='displayexp', label='Take exponential transformation', value=FALSE),
                            h6("Reference"),
                            checkboxInput(inputId='displayref', label='Display no subgroup effect outcome',
                                          value=TRUE)
                            )
                 ),
                 h6("Select subgroups to display in analysis results"),
                 DT::dataTableOutput('subgrpsel')
             ))
}

##download report
panel.download <- function() {
    rst <- list(wellPanel(h4("Summary"),
                   htmlOutput("rstsummary"),
                   DT::dataTableOutput('rsttbl')));

    has.rd <- requireNamespace("rmarkdown", quietly = TRUE);
    has.pd <- requireNamespace("pander", quietly = TRUE);

    if (has.rd & has.pd) {
        rst <- c(rst,
                 list(wellPanel(h4('Download the analysis report'),
                           radioButtons('format', '', c('PDF', 'HTML', 'Word')),
                           downloadButton('btnDload')))
                 );
    } else {
        tmp.msg <- NULL;
        if (!has.rd)
            tmp.msg <- "rmarkdown";
        if (!has.pd)
            tmp.msg <- c(tmp.msg,  "pander");
        msg <- paste('R package ',
                     paste(tmp.msg, collapse=" and "),
                     ' needed for generating the analysis report. Please install.');
        rst <- c(rst,
                 list(wellPanel(msg.box(msg, "warning")))
                 );
    }
    rst
}

##define the main tabset for beans
tab.main <- function() {
    panels <- list(position = "above", type = "pills",
                   id="mainpanel",
                   tab.start()
                  ,tab.upload()
                  ,tabPanel("Subgroup Specification", uiOutput("uiSubgroup"))
                  ,tab.mcmc()
                  ,tabPanel("Bayesian Analysis", uiOutput("uiAnalysis"))
                  ,tabPanel("Toolbox", uiOutput("uiTool"))
                  ,tabPanel("Report", uiOutput("uiReport"))
                   );

    ##if (1 == userLog$version) {
    ##    panels <- c(panels, list(tab.mcmc())); ##tab.model(),
    ##}

    ##generate
    do.call(tabsetPanel, panels);
}


##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------

##check subgroup specification for raw data
chk.raw.subgrp <- function(dset, resp.type, trt, resp, censor, cov, levels) {

    msg <- NULL;

    if (trt == resp)
        msg <- paste(msg, "<li>Column used multiple times </li>");

    if (any(is.na(dset[,c(trt, resp, censor, cov)])))
        msg <- paste(msg, "<li>Missing value found </li>");

    if (is.null(cov)) {
        msg <- paste(msg, "<li>No covariates selected</li>");
    }

    if (!is.null(levels)) {
        if (any(1 >= levels))
            msg <- paste(msg, "<li>Some covariate has less than two levels</li>");
    };

    if (length(unique(dset[, trt])) != 2)
        msg <- paste(msg, "<li>There are not two arms of treatment </li>");

    if (resp.type == RESP.TYPE[3]) {
        if (trt == censor | resp == censor)
            msg <- paste(msg, "<li>Column used multiple times </li>");

        if (length(unique(dset[, censor])) != 2)
            msg <- paste(msg, "<li>Censoring indicator has invalid value </li>");
    }

    if (resp.type == RESP.TYPE[2]) {
        if (length(unique(dset[, resp])) != 2)
            msg <- paste(msg, "<li>Response is not binary </li>");
    }



    ##if (nrow(unique(dset[, c(trt, cov)])) != nrow(dset))
    ##    msg <- paste(msg, "<li>Subgroups are not unique</li>");

    ##return
    if (is.null(msg)) {
        rst <- NULL;
    } else {
        rst <- "Subgroup specification is not valid. Please check the following:";
        rst <- paste(rst, "<ul>", msg, "</ul>");
    }

    rst
}

##check subgroup specification for summary trt effect data
chk.summary.subgrp <- function(dset, est, var, cov, levels) {

    ##print(unique(dset[, cov]));
    msg <- NULL;
    if (any(is.na(dset[,c(est, var, cov)])))
        msg <- paste(msg, "<li>Missing value found </li>");

    if (is.factor(dset[,var])|is.factor(dset[,est]))
        msg <- paste(msg, "<li>Estimate and variance columns need to be numerical values </li>");

    if (any(as.numeric(dset[,var]) <= 0))
        msg <- paste(msg, "<li>Negative variance value found </li>");

    if (is.null(cov)) {
        msg <- paste(msg, "<li>No covariates selected</li>");
    };

    if (!is.null(levels)) {
        if (any(1 >= levels))
            msg <- paste(msg, "<li>Some covariate has less than two levels</li>");
    };

    ##else if (nrow(unique(dset[, cov, drop=FALSE])) != nrow(dset)) {
    ##    msg <- paste(msg, "<li>Subgroups are not unique</li>");
    ##}

    ##return
    if (is.null(msg)) {
        rst <- NULL;
    } else {
        rst <- "Subgroup specification is not valid. Please check the following:";
        rst <- paste(rst, "<ul>", msg, "</ul>");
    }

    rst
}


observe ({
    inFile <- input$userdata;
    if (!is.null(inFile)){
        userLog$data <- read.csv(inFile$datapath,
                                    header=input$header,
                                    sep=input$sep,
                                    quote=input$quote)
    }
})

observe ({
    if (is.null(input$btnExample))
        return(NULL);

    if (0 == input$btnExample)
        return(NULL);

    userLog$data <- solvd.sub;
})

##get subgroup covariates
get.sub.cov <- reactive({
    ##if (-1 == userLog$SubgroupValid)
    ##    return(NULL);

    if (DATA.FORMAT[1] == input$dataformat) {
        rst <- input$selsubcov;
    } else {
        rst <- input$selcov;
    }

    rst
})

##get levels of covariates
get.sub.cov.levels <- reactive({
    c.covs <- get.sub.cov();
    c.data <- get.data();

    if (is.null(c.covs) | is.null(c.data))
        return(NULL);

    rst <- NULL;
    for (i in 1:length(c.covs)) {
        tmp.l <- length(unique(c.data[, c.covs[i]]));
        rst   <- c(rst, tmp.l);
    }

    rst
})

##read data from user input
get.data <- reactive({
    userLog$data;
})

##get subgroup effects
get.subgrp.raw <- function() {
    c.data <- get.data();
    if (is.null(c.data)        |
        is.null(input$selcov)  |
        is.null(input$selresp) |
        is.null(input$seltrt)  |
        is.null(input$resptype)
        )
        return(NULL);

    resptype <- c("continuous", "binary", "survival")[which(RESP.TYPE == input$resptype)];
    rst <- r.get.subgrp.raw(c.data, input$selresp, input$seltrt,
                            input$selcov, input$selcensor, resptype);
}

##get sub group effects from subgroup summary file
get.subgrp.sub <- function() {
    c.data <- get.data();
    if (is.null(c.data) |
        is.null(input$selsubcov) |
        is.null(input$selsube)|
        is.null(input$selsubvar)
        )
        return();

    rst <- r.subgrp.effect(c.data, input$selsube, input$selsubvar, input$selsubcov);
}

##get subgroup data
get.subgrp <- reactive({
    if (is.null(input$btnSub))
        return(NULL);

    if (0 == input$btnSub)
        return(NULL);

    ##if (-1 == userLog$SubgroupValid)
    ##    return(NULL);
    isolate({
        rst <- tryCatch(
            {
                if (input$dataformat == DATA.FORMAT[2]) {
                    ##Create a Progress object
                    progress <- shiny::Progress$new(session, min=0, max=1);
                    progress$set(message = "Get subgroup treatment effects...", value=1);
                    on.exit(progress$close());
                    rst <- get.subgrp.raw();
                } else {
                    rst <- get.subgrp.sub();
                }
            },
            error=function(cond) {
                print(cond);
                return(NULL);
            }
        );

        rst});
})

##get selected subgroups
get.subgrp.selected <- reactive({

    input$btnAna;
    input$subgrpsel_rows_selected;

    isolate({
        dat.sub <- get.subgrp();
        if (is.null(dat.sub))
            return(NULL);

        if (is.null(input$subgrpsel_rows_selected)) {
            rst <- 1:nrow(dat.sub);
        } else {
            rst <- sort(input$subgrpsel_rows_selected);
        }
    })
    rst
})

##get model selected for analysis
##no subgroup effect model is always selected
get.ana.models <- reactive({

    ## if (0 == userLog$version) {
    ##     rst <- ALL.MODELS;
    ## } else {
    ##     rst <- input$selmdl;
    ## }

    rst <- ALL.MODELS[1];
    for (i in 2:length(ALL.MODELS)) {
        cur.in <- paste("chkM", i, sep="");
        if (input[[cur.in]]) {
            rst <- c(rst, ALL.MODELS[i]);
        }
    }

    rst
})

##get mcmc setup
get.mcmc.par <- reactive({
    rst <- list(iter=input$mcmciter,
                warmup=input$mcmcburnin,
                thin=input$mcmcthin,
                seed=input$mcmcseed,
                algorithm=input$mcmcalg);
    rst
})

##get prior parameters
get.par.prior <- reactive({
    rst <- list(mdl1=list(vtau=input$m1vtau),
                mdl2=list(vtau=input$m2vtau),
                mdl3=list(vtau=input$m3vtau, vgamma=input$m3vgamma),
                mdl4=list(vtau=input$m4vtau, vw=input$m4vw),
                mdl5=list(vtau=input$m5vtau, vw=input$m5vw, vgamma=input$m5vgamma),
                mdl6=list(vtau=input$m6vtau, vw=input$m6vw),
                mdl7=list(vtau=input$m7vtau, vw=input$m7vw));

    rst <- mapply(function(x) {c(x, list(vrange=input$rangelogvar))},
                  rst, SIMPLIFY=FALSE);
    rst
})

##trace plot height
get.traceplot.height <- function(smps, n.eachrow=2, row.height=100) {
    n   <- dim(smps)[3];
    rst <- n/n.eachrow*row.height;
}


##get dic for all models
get.all.dic <- reactive({
    ##model results
    mrst <- ana.rst();

    if (is.null(mrst))
        return(NULL);

    rst <- NULL;
    for (i in 1:length(mrst)) {
        rst <- c(rst, mrst[[i]]$dic);
    }

    rst
})


##get anoint analysis results
get.anoint <- reactive({

    c.data <- get.data();

    if (is.null(c.data)        |
        is.null(input$selcov)  |
        is.null(input$selresp) |
        is.null(input$seltrt)  |
        is.null(input$resptype)
        )
        return(NULL);

    print(input$selcensor);

    rst <- r.get.anoint(c.data, input$selresp, input$seltrt,
                        input$selcov, input$selcensor, input$resptype);
})



##apply anoint to subject level data
r.get.anoint <- function(data.all, var.resp, var.trt, var.cov, var.censor, resptype) {

    if (resptype == RESP.TYPE[3]) {
        ##survival
        resp       <- paste("survival::Surv(", var.resp, ",", var.censor, ")", sep="");
        mdl.family <- "coxph";
    } else {
        ##continuous or binary
        resp <- var.resp;
        if (resptype == RESP.TYPE[1]) {
            mdl.family <- "gaussian";
        } else if (resptype == RESP.TYPE[2]) {
            mdl.family <- "binomial";
        }
    }

    ## interaction models
    mdl.formula <- as.formula(paste(resp, "~ (",
                                    paste(var.cov, collapse="+"), ")*",
                                    var.trt, sep=""));
    ano.rst <- do.call(anoint::anoint,
                       list(formula=mdl.formula,
                            data=data.all,
                            family=mdl.family));

    ## proportional interaction models
    mdl.fml.2 <- as.formula(paste(resp, "~ ", paste(var.cov, collapse="+"), sep=""));
    pim.rst   <- do.call(anoint::pim.fit,
                         list(formula=mdl.fml.2,
                              data=data.all,
                              family=mdl.family, trt=var.trt));

    ##return
    list(obo=anoint::obo(ano.rst),
         uim=anoint::uim(ano.rst),
         pim=pim.rst);
}
