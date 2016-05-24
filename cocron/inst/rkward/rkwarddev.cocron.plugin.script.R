require("rkwarddev")

local({
  # set the output directory to overwrite the actual plugin
  output.dir <- "~/cocron/" # tempdir()
  overwrite <- TRUE

  about.node <- rk.XML.about(
    name="cocron",
    author=person(given="Birk", family="Diedenhofen", email="mail@birkdiedenhofen.de", role=c("aut","cre")),
    about=list(
      desc="Compare two or more alpha coefficients based on either dependent or independent groups of individuals",
      version="1.0-1",
      date="2016-03-11", # Sys.Date(),
      url="http://comparingcronbachalphas.org",
      license="GPL"
    )
  )

  dependencies.node <- rk.XML.dependencies(
    dependencies=list(
      rkward.min="0.6.1",
      R.min="2.15"
    )
  )

  data.input <- rk.XML.radio("Data input",
    id.name="data_input",
    options=list(
      "Calculate and compare Cronbach alphas from raw data"=c(val="raw.data"),
      "Enter alpha coefficients manually"=c(val="manual")
    )
  )

  groups <- rk.XML.radio("The Cronbach alphas are based on",
    id.name="groups",
    options=list(
      "two independent groups"=c(val="indep"),
      "two dependent groups (e.g., same group)"=c(val="dep")
    )
  )

  alpha.count <- rk.XML.spinbox("Number of alpha coefficients to be compared:", min=2, max=7, id.name="alpha_count", initial=2, real=FALSE)

  manual.and.indep.groups.alpha <- rk.XML.matrix(" ", mode="real", columns=1, min=0, max=1, allow_user_resize_rows=FALSE, allow_user_resize_columns=FALSE, fixed_width=TRUE, horiz_headers="Alpha", id.name="manual_and_indep_groups_alpha")
  manual.and.indep.groups.n <- rk.XML.matrix(" ", mode="integer", columns=1, min=0, allow_user_resize_rows=FALSE, allow_user_resize_columns=FALSE, fixed_width=TRUE, horiz_headers="Sample size", id.name="manual_and_indep_groups_n")
  manual.and.indep.groups.i <-rk.XML.matrix(" ", mode="integer", columns=1, min=0, allow_user_resize_rows=FALSE, allow_user_resize_columns=FALSE, fixed_width=TRUE, horiz_headers="Item count", id.name="manual_and_indep_groups_i")

  manual.and.dep.groups.alpha <- rk.XML.matrix(" ", mode="real", columns=1, min=0, max=1, allow_user_resize_rows=FALSE, allow_user_resize_columns=FALSE, fixed_width=TRUE, horiz_headers="Alpha", id.name="manual_and_dep_groups_alpha")
  manual.and.dep.groups.n <- rk.XML.spinbox("Sample size:", min=0, id.name="manual_and_dep_groups_n", initial=100, real=FALSE)
  manual.and.dep.groups.i <-rk.XML.spinbox("Item count:", min=0, id.name="manual_and_dep_groups_i", initial=10, real=FALSE)
  manual.and.dep.groups.r <- rk.XML.matrix(" ", mode="real", min=-1, max=1, allow_user_resize_rows=FALSE, allow_user_resize_columns=FALSE, id.name="manual_and_dep_groups_r")

  var.sel <- rk.XML.varselector("Data", id.name="vars")
  raw.data <- rk.XML.varslot("Select at least two data.frames/matrices", source=var.sel, classes=c("data.frame", "matrix"), multi=TRUE, min=2, id.name="raw_data_alpha")

  standardized.alpha <- rk.XML.cbox("Calculate standardized Cronbach's alpha", value = "true", un.value = "false", chk = FALSE, id.name = "standardized_alpha")
  standardized.alpha.frame <- rk.XML.frame(standardized.alpha, label="Should the scores be standardized before calculating the alpha coefficients?")

  los <- rk.XML.spinbox("Level of significance:", min=0, max=1, id.name="los", initial=.05)
  los.frame <- rk.XML.frame(los, label="Please choose a level of significance for the test of significance:")

  conf.level <- rk.XML.spinbox("Level of confidence", min=0, max=1, id.name="conf_int", initial=.95)
  conf.level.frame <- rk.XML.frame(conf.level, label="Please choose a confidence level for the confidence intervals:")

  wizard.raw.data.input.page <- rk.XML.page(
    rk.XML.row(
      rk.XML.text("Please provide the raw data the Cronbach alphas should be calculated from:<br />", id.name="please_provide_the_raw_data"),
      id.name="row_instruction_to_provide_raw_data"
    ),
    rk.XML.row(
      rk.XML.col(var.sel),
      rk.XML.col(raw.data),
      id.name="row_input_raw_data"
    )
  )

  wizard.manual.data.input.page <- rk.XML.page(
    rk.XML.row(
      rk.XML.text("Please provide the Cronbach alphas you want to compare and the sample sizes and number of items they are based on:<br />", id.name="please_provide_the_cronbach_alphas"),
      id.name="row_instruction_to_input_cronbach_alpha_and_sample_size"
    ),
    rk.XML.row(
      alpha.count,
      id.name="row_input_cronbach_alpha_count"
    ),
    rk.XML.row(
      rk.XML.col(manual.and.dep.groups.n),
      rk.XML.col(manual.and.dep.groups.i),
      id.name="row_input_dep_groups_n_and_i"
    ),
    rk.XML.row(
      rk.XML.row(
        manual.and.indep.groups.alpha,
        id.name="row_input_indep_groups_alphas"
      ),
      rk.XML.row(
        rk.XML.col(manual.and.indep.groups.n),
        rk.XML.col(manual.and.indep.groups.i),
        id.name="row_input_indep_groups_n_and_i"
      ),
      rk.XML.row(
        manual.and.dep.groups.alpha,
        id.name="row_input_dep_groups_alphas"
      ),
      rk.XML.stretch()
    )
  )

  wizard.correlations.page <- rk.XML.page(
    rk.XML.text("Please provide the correlations between the underlying scores of the Cronbach alphas:<br />"),
    manual.and.dep.groups.r
  )

  wizard <- rk.XML.wizard(
    label="Comparing Cronbach alphas",
    rk.XML.page(
      rk.XML.text("Are the Cronbach alphas based on independent or on dependent groups? (If the data were taken from measurements of the same individuals, the groups are dependent.)<br />"),
      groups,
      rk.XML.stretch()
    ),
    rk.XML.page(
      rk.XML.text("Do you want to calculate and compare alpha coefficients from raw data or do you want to type the alpha coefficients in manually?"),
      data.input,
      rk.XML.stretch()
    ),
    wizard.raw.data.input.page,
    wizard.manual.data.input.page,
    wizard.correlations.page,
    rk.XML.page(
      standardized.alpha.frame,
      los.frame,
      conf.level.frame,
      rk.XML.stretch()
    )
  )

  logic <- rk.XML.logic(
    ## convert single
    raw.data.input.convert <- rk.XML.convert(sources=list(string=data.input), mode=c(equals="raw.data"), id.name="raw_data_input_convert"),
    manual.input.convert <- rk.XML.convert(sources=list(string=data.input), mode=c(equals="manual"), id.name="manual_input_convert"),

    indep.groups.convert <- rk.XML.convert(sources=list(string=groups), mode=c(equals="indep"), id.name="indep_groups_convert"),
    dep.groups.convert <- rk.XML.convert(sources=list(string=groups), mode=c(equals="dep"), id.name="dep_groups_convert"),

    ## convert multiple
    raw.data.and.indep.groups.convert <- rk.XML.convert(sources=list(raw.data.input.convert, indep.groups.convert), mode=c(and=""), id.name="raw_data_and_indep_groups_convert"),
    raw.data.and.dep.groups.convert <- rk.XML.convert(sources=list(raw.data.input.convert, dep.groups.convert), mode=c(and=""), id.name="raw_data_and_dep_groups_convert"),

    manual.and.indep.groups.convert <- rk.XML.convert(sources=list(manual.input.convert, indep.groups.convert), mode=c(and=""), id.name="manual_and_indep_groups_convert"),
    manual.and.dep.groups.convert <- rk.XML.convert(sources=list(manual.input.convert, dep.groups.convert), mode=c(and=""), id.name="manual_and_dep_groups_convert"),

    ## connect
    rk.XML.connect(governor=raw.data.input.convert, client=wizard.raw.data.input.page, set="visible"),
    rk.XML.connect(governor=manual.input.convert, client=wizard.manual.data.input.page, set="visible"),
    rk.XML.connect(governor=manual.and.dep.groups.convert, client=wizard.correlations.page, set="visible"),

    rk.XML.connect(governor=manual.and.indep.groups.convert, client=manual.and.indep.groups.alpha, set="visible"),
    rk.XML.connect(governor=manual.and.indep.groups.convert, client=manual.and.indep.groups.n, set="visible"),
    rk.XML.connect(governor=manual.and.indep.groups.convert, client=manual.and.indep.groups.i, set="visible"),

    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.alpha, set="visible"),
    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.n, set="visible"),
    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.i, set="visible"),

    rk.XML.connect(governor=raw.data.input.convert, client=standardized.alpha.frame, set="visible"),

    rk.XML.connect(governor="alpha_count.int", client=manual.and.dep.groups.r, set="rows"),
    rk.XML.connect(governor="alpha_count.int", client=manual.and.dep.groups.r, set="columns"),

    rk.XML.connect(governor="alpha_count.int", client=manual.and.indep.groups.alpha, set="rows"),
    rk.XML.connect(governor="alpha_count.int", client=manual.and.indep.groups.n, set="rows"),
    rk.XML.connect(governor="alpha_count.int", client=manual.and.indep.groups.i, set="rows"),

    rk.XML.connect(governor="alpha_count.int", client=manual.and.dep.groups.alpha, set="rows"),

    ## require
    rk.XML.connect(governor=raw.data.input.convert, client=raw.data, set="required"),

    rk.XML.connect(governor=manual.and.indep.groups.convert, client=manual.and.indep.groups.alpha, set="required"),
    rk.XML.connect(governor=manual.and.indep.groups.convert, client=manual.and.indep.groups.n, set="required"),
    rk.XML.connect(governor=manual.and.indep.groups.convert, client=manual.and.indep.groups.i, set="required"),

    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.alpha, set="required"),
    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.n, set="required"),
    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.i, set="required"),
    rk.XML.connect(governor=manual.and.dep.groups.convert, client=manual.and.dep.groups.r, set="required")
  )

  JS.calc <- rk.paste.JS(
    ## raw data input
    ite(id(data.input, " == 'raw.data'"), echo("result <- cocron(")),

    raw.data.vars <- rk.JS.vars(raw.data, join = ", "),
    ite(id(data.input, " == 'raw.data' && ", groups, " == 'indep'"), echo("data=list(", raw.data.vars, "), dep=FALSE")),
    ite(id(data.input, " == 'raw.data' && ", groups, " == 'dep'"), echo("data=list(", raw.data.vars, "), dep=TRUE")),

    ite(id(data.input, " == 'raw.data' && ", standardized.alpha, " == 'true'"), echo(", standardized=TRUE")),

    ## manual input
    ite(id(data.input, " == 'manual'"), echo("result <- cocron.n.coefficients(")),

    ite(id(data.input, " == 'manual' && ", groups, " == 'indep'"), echo("alpha=as.vector(", manual.and.indep.groups.alpha, "), n=as.vector(", manual.and.indep.groups.n, "), items=as.vector(", manual.and.indep.groups.i, "), dep=FALSE")),
    ite(id(data.input, " == 'manual' && ", groups, " == 'dep'"), echo("alpha=as.vector(", manual.and.dep.groups.alpha, "), n=", manual.and.dep.groups.n, ", items=", manual.and.dep.groups.i, ", dep=TRUE, r=", manual.and.dep.groups.r)),

    echo(", los=",  los, ", conf.level=", conf.level, ")\n"),
    level=2
  )

  JS.print <- rk.paste.JS(echo("rk.print(result)\n"), level=2)

  rkh <- list(
    summary=rk.rkh.summary(text="Compare two or more alpha coefficients based on either dependent or independent groups of individuals."),
    sections=list(
      groups=rk.rkh.section(title="Step 1", text="Indicate if the alpha coefficients are based on independent or dependent groups."),
      data.input=rk.rkh.section(title="Step 2", text="Choose whether you want to calculate the alpha coefficients from raw data or type them in manually."),
      data=rk.rkh.section(title="Step 3", text="Enter the data."),
      test.hypothesis=rk.rkh.section(title="Step 4", text="If the alpha coefficients are calculated from raw data, choose whether the item scores should be standardized. Select a level of significance for the significance test and a level of confidence for the confidence intervals."),
      code=rk.rkh.section(title="Step 5", text="Copy and paste the generated code to your R script or directly run the code.")
    )
  )

  plugin.dir <<- rk.plugin.skeleton(
    about.node,
    path=output.dir,
    provides=c("logic", "wizard"),
    xml=list(logic=logic, wizard=wizard),
    rkh=rkh,
    js=list(
      require="cocron",
      results.header="Comparing Cronbach alphas",
      calculate=JS.calc,
      printout=JS.print
    ),
    pluginmap=list(name="Comparing Cronbach alphas", hierarchy=list("analysis", "Classical test theory")),
    dependencies=dependencies.node,
    guess.getter=TRUE,
    load=TRUE,
#     edit=TRUE,
#     show=TRUE,
    overwrite=overwrite
  )
})
