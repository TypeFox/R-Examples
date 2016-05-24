######################################################################################################################

# Function: CreateReportStructure.
# Argument: Results returned by the CSE function and tables produced by the CreateTableStructure function.
# Description: This function is used to create a report

CreateReportStructure = function(evaluation, presentation.model){

  # Number of scenario
  n.scenario = nrow(evaluation$analysis.scenario.grid)

  # Number of design parameter
  n.design = max(evaluation$analysis.scenario.grid$design.parameter)

  # Number of outcome parameter
  n.outcome = max(evaluation$analysis.scenario.grid$outcome.parameter)

  # Number of sample size or event set
  n.sample.size = max(evaluation$analysis.scenario.grid$sample.size)

  # Number of multiplicity adjustment
  n.multiplicity.adjustment = max(evaluation$analysis.scenario.grid$multiplicity.adjustment)

  # Empty report object
  report = list()

  # Empty section list
  section = list()

  # Empty subsection list
  subsection = list()

  # Empty subsubsection list
  subsubsection = list()

  # Empty subsubsubsection list
  subusbsubsection = list()

  # Get the label
  custom.label = presentation.model$custom.label

  # Sample size label
  custom.label.sample.size = list()
  if (any(unlist(lapply(custom.label, function(x) (x$param %in% c("sample.size","event")))))) {
    custom.label.sample.size$label = custom.label[[which(unlist(lapply(custom.label, function(x) (x$param %in% c("sample.size","event")))))]]$label
    custom.label.sample.size$custom = TRUE

  } else {
    if (any(!is.na(evaluation$data.structure$sample.size.set))) custom.label.sample.size$label = paste("Sample size", 1:n.sample.size)
    else if (any(!is.na(evaluation$data.structure$event.set))) custom.label.sample.size$label = paste("Event", 1:n.sample.size)
    custom.label.sample.size$custom = FALSE
  }

  # Outcome parameter label
  custom.label.outcome.parameter = list()
  if (any(unlist(lapply(custom.label, function(x) (x$param == "outcome.parameter"))))) {
    custom.label.outcome.parameter$label = custom.label[[which(unlist(lapply(custom.label, function(x) (x$param == "outcome.parameter"))))]]$label
    custom.label.outcome.parameter$custom = TRUE
  } else {
    custom.label.outcome.parameter$label = paste("Outcome", 1:n.outcome)
    custom.label.outcome.parameter$custom = FALSE
  }

  # Multiplicity adjustment label
  custom.label.multiplicity.adjustment = list()
  if (any(unlist(lapply(custom.label, function(x) (x$param == "multiplicity.adjustment"))))) {
    custom.label.multiplicity.adjustment$label = custom.label[[which(unlist(lapply(custom.label, function(x) (x$param == "multiplicity.adjustment"))))]]$label
    custom.label.multiplicity.adjustment$custom = TRUE
  } else {
    custom.label.multiplicity.adjustment$label = paste("Multiplicity adjustment scenario", 1:n.multiplicity.adjustment)
    custom.label.multiplicity.adjustment$custom = FALSE
  }

  # Design parameter label
  custom.label.design.parameter = list()
  if (any(unlist(lapply(custom.label, function(x) (x$param == "design.parameter"))))) {
    custom.label.design.parameter$label = custom.label[[which(unlist(lapply(custom.label, function(x) (x$param == "design.parameter"))))]]$label
    custom.label.design.parameter$custom = TRUE
  } else {
    custom.label.design.parameter$label = paste("Design", 1:n.design)
    custom.label.design.parameter$custom = FALSE
  }

  # Create a summary table for the design
  if (!is.null(evaluation$data.structure$design.parameter.set)) table.design = CreateTableDesign(evaluation$data.structure, custom.label.design.parameter$label)

  # Create a summary table for the sample size
  table.sample.size = CreateTableSampleSize(evaluation$data.structure, custom.label.sample.size$label)

  # Create a summary table for the outcome parameters
  outcome.information = CreateTableOutcome(evaluation$data.structure, custom.label.outcome.parameter$label)
  outcome.dist.name = outcome.information[[1]]
  table.outcome =  outcome.information[[2]]

  # Create a summary table for the tests
  table.test = CreateTableTest(evaluation$analysis.structure)

  # Create a summary table for the statistics
  if (!is.null(evaluation$analysis.structure$statistic)) table.statistic = CreateTableStatistic(evaluation$analysis.structure)

  # Create a summary table for the results, according to the section/subsection requested by the user
  result.structure = CreateTableStructure(evaluation, presentation.model, custom.label.sample.size, custom.label.design.parameter, custom.label.outcome.parameter, custom.label.multiplicity.adjustment)

  # Get information on the multiplicity adjustment
  mult.adj.desc = list()
  if (!is.null(evaluation$analysis.structure$mult.adjust)){
    for (mult in 1:n.multiplicity.adjustment) {
      mult.adjust.temp = list()
      # Number of multiplicity adjustment within each mult.adj scenario
      n.mult.adj.sc = length(evaluation$analysis.structure$mult.adjust[[mult]])
      for (j in 1:n.mult.adj.sc){
        if (!is.na(evaluation$analysis.structure$mult.adjust[[mult]][[j]]$proc)){
          dummy.function.call = list("Description", evaluation$analysis.structure$mult.adjust[[mult]][[j]]$par, unlist(evaluation$analysis.structure$mult.adjust[[mult]][[j]]$tests))
          analysis.mult.desc = do.call(evaluation$analysis.structure$mult.adjust[[mult]][[j]]$proc, list(rep(0,length(unlist(evaluation$analysis.structure$mult.adjust[[mult]][[j]]$tests))),dummy.function.call))
          mult.adjust.temp[[j]] = list(desc = analysis.mult.desc[[1]], tests = paste0("{",paste0(unlist(evaluation$analysis.structure$mult.adjust[[mult]][[j]]$tests),collapse=", "),"}"),par = analysis.mult.desc[[2]])
        }
        else {
          mult.adjust.temp[[j]] = list(desc = "No adjustment", tests=NULL, par=NULL)
        }
      }
      mult.adj.desc[[mult]] = mult.adjust.temp
    }
  } else {
    mult.adj.desc = NA
  }

  # Section 1: General information
  ##################################

  # Items included in Section 1, Subsection 1
  # Item's type is text by default
  item1 = list(label = "", value = paste0("This report was generated by ", presentation.model$project$username, " using the Mediana package. For more information about the Mediana package, see http://gpaux.github.io/Mediana."))
  item2 = list(label = "Project title:", value = presentation.model$project$title)
  item3 = list(label = "Description:", value = presentation.model$project$description)
  item4 = list(label = "Random seed:", value = evaluation$sim.parameters$seed)
  item5 = list(label = "Number of simulations:", value = evaluation$sim.parameters$n.sims)
  item6 = list(label = "Number of cores:", value = evaluation$sim.parameters$proc.load)
  item7 = list(label = "Start time:", value = evaluation$timestamp$start.time)
  item8 = list(label = "End time:", value = evaluation$timestamp$end.time)
  item9 = list(label = "Duration (mins):", value = format(round(evaluation$timestamp$duration, digits = 2), digits = 2, nsmall = 2))

  # Create a subsection (set the title to NA to suppress the title)
  subsection[[1]] = list(title = "Project information", item = list(item1, item2, item3))
  # Create a subsection (set the title to NA to suppress the title)
  subsection[[2]] = list(title = "Simulation parameters", item = list(item4, item5, item6, item7, item8, item9))

  # Create the header section (set the title to NA to suppress the title)
  section[[1]] = list(title = "General information", subsection = subsection)

  # Section 2: Data model #
  #########################
  n.subsection = 0

  # Empty subsection list
  subsection = list()

  # Empty subsubsection list
  subsubsection = list()

  # Empty subsubsubsection list
  subusbsubsection = list()

  #Design parameters
  if (!is.null(evaluation$data.structure$design.parameter.set)){
    n.subsection = n.subsection + 1
    item1 = list(label = "Number of design parameter sets: ",
                 value = n.design
    )
    item2 = list(label = "Design",
                 value =  table.design[,2:length(table.design)],
                 param = list(groupedheader.row = list(values = c("", "Enrollment", "", "", "Dropout"), colspan = c(1, 3, 1, 1, 2))),
                 type = "table"
    )

    # Create a subsection
    subsection[[n.subsection]] = list(title = "Design", item = list(item1, item2))
  }

  #Sample size
  if (any(!is.na(evaluation$data.structure$sample.size.set))) {
    item1 = list(label = "Number of samples:",
                 value = length(evaluation$data.structure$id),
                 type = "text"
    )
    item2 = list(label = "Number of sample size sets:",
                 value = n.sample.size,
                 type = "text"
    )
    item3 = list(label = "Sample size",
                 value = table.sample.size[,2:ncol(table.sample.size)],
                 param = list(span.columns = "Sample size set"),
                 type = "table"
    )
    # Create a subsection
    subsection[[n.subsection+1]] = list(title = "Sample size", item = list(item1, item2, item3))
  }
  #Event
  if (any(!is.na(evaluation$data.structure$event.set))) {
    item1 = list(label = "Number of samples:",
                 value = length(evaluation$data.structure$id),
                 type = "text"
    )
    item2 = list(label = "Randomization ratio:",
                 value = paste0("(",paste0(evaluation$data.structure$rando.ratio, collapse = ":"),")"),
                 type = "text"
    )
    item3 = list(label = "Number of event sets:",
                 value = n.sample.size,
                 type = "text"
    )
    item4 = list(label = "Event",
                 value = table.sample.size[,2:ncol(table.sample.size)],
                 type = "table"
    )
    # Create a subsection
    subsection[[n.subsection+1]] = list(title = "Number of events", item = list(item1, item2, item3, item4))
  }

  # Outcome distribution
  item1 = list(label = "Number of outcome parameter sets:",
               value = n.outcome,
               type = "text"
  )
  item2 = list(label = "Outcome distribution:",
               value = outcome.dist.name,
               type = "text"
  )
  item3 = list(label = "Outcome parameter",
               value = table.outcome[,2:length(table.outcome)],
               param = list(span.columns = "Outcome parameter set"),
               type = "table"
  )
  # Create a subsection
  subsection[[n.subsection+2]] = list(title = "Outcome distribution", item = list(item1, item2, item3))

  section[[2]] = list(title = "Data model", subsection = subsection)

  # Section 3: Analysis model
  ###########################
  n.subsection = 0

  # Empty subsection list
  subsection = list()

  # Empty subsection list
  subsection = list()

  # Empty subsubsection list
  subsubsection = list()

  # Empty subsubsubsection list
  subusbsubsection = list()

  # Test
  if (!is.null(evaluation$analysis.structure$test)){
    n.subsection = n.subsection + 1
    item1 = list(label = "Number of tests/null hypotheses: ",
                 value = length(evaluation$analysis.structure$test)
    )
    item2 = list(label = "Tests",
                 value =  table.test,
                 type = "table"
    )

    # Create a subsection
    subsection[[n.subsection]] = list(title = "Tests", item = list(item1, item2))

  }

  # Statistic
  if (!is.null(evaluation$analysis.structure$statistic)){
    n.subsection = n.subsection + 1
    item1 = list(label = "Number of descriptive statistics: ",
                 value = length(evaluation$analysis.structure$statistic)
    )
    item2 = list(label = "Statistics",
                 value =  table.statistic,
                 type = "table"
    )

    # Create a subsection
    subsection[[n.subsection]] = list(title = "Statistics", item = list(item1, item2))

  }

  # Multiplicity adjustment
  if (!is.null(evaluation$analysis.structure$mult.adjust)){
    n.subsection = n.subsection + 1
    subsubsection = list()
    for (mult in 1:n.multiplicity.adjustment) {
      # Number of multiplicity adjustment within each mult.adj scenario
      n.mult.adj.sc = length(mult.adj.desc[[mult]])
      subsubsubsection = list()
      for (j in 1:n.mult.adj.sc){
        item = list()
        ind.item = 1
        item[[ind.item]] = list(label = "Procedure:",
                                value = mult.adj.desc[[mult]][[j]]$desc[[1]]
        )
        if (!is.null(mult.adj.desc[[mult]][[j]]$tests)){
          ind.item = ind.item + 1
          item[[ind.item]] = list(label = "Tests:",
                                  value =  mult.adj.desc[[mult]][[j]]$tests
          )
        }
        if (!is.null(mult.adj.desc[[mult]][[j]]$par)){
          ind.item = ind.item + 1
          if (length(mult.adj.desc[[mult]][[j]]$par)>1) {
            item[[ind.item]] = list(label = "Parameters:",
                                    value =  ""
            )
            for (k in 1:length(mult.adj.desc[[mult]][[j]]$par)){
              ind.item = ind.item + 1
              if (!is.data.frame(mult.adj.desc[[mult]][[j]]$par[[k]])) {
                item[[ind.item]] = list(label = "",
                                        value = mult.adj.desc[[mult]][[j]]$par[[k]],
                                        type = "text"
                )
              }
              else if (is.data.frame(mult.adj.desc[[mult]][[j]]$par[[k]])) {
                item[[ind.item]] = list(label = "Parameters",
                                        value = mult.adj.desc[[mult]][[j]]$par[[k]],
                                        type = "table"
                )
              }
            }
          }
          else {
            if (!is.data.frame(mult.adj.desc[[mult]][[j]]$par[[1]])) {
              item[[ind.item]] = list(label = "Parameters:",
                                      value = mult.adj.desc[[mult]][[j]]$par[[1]],
                                      type = "text"
              )
            }
            else if (is.data.frame(mult.adj.desc[[mult]][[j]]$par[[1]])) {
              item[[ind.item]] = list(label = "Parameters:",
                                      value = mult.adj.desc[[mult]][[j]]$par[[1]],
                                      type = "table"
              )
            }
          }
        }
        if (n.mult.adj.sc>1) {
          subsubsubsection[[j]] = list(title = paste0("Multiplicity adjustment procedure ",j), item = item)
        }
      }
      if (n.mult.adj.sc>1) {
        subsubsection[[mult]] = list(title = custom.label.multiplicity.adjustment$label[mult], subsubsubsection = subsubsubsection)
      }
      else if (!is.null(evaluation$analysis.structure$mult.adjust) & n.multiplicity.adjustment>1){
        subsubsection[[mult]] = list(title = custom.label.multiplicity.adjustment$label[mult], item = item)
      }
    }
    if (n.mult.adj.sc>1) {
      subsection[[n.subsection]] = list(title = "Multiplicity adjustment", subsubsection = subsubsection )
    }
    else if (!is.null(evaluation$analysis.structure$mult.adjust) & n.multiplicity.adjustment>1){
      subsection[[n.subsection]] = list(title = "Multiplicity adjustment", subsubsection = subsubsection )
    }
    else if (!is.null(evaluation$analysis.structure$mult.adjust) & n.multiplicity.adjustment==1){
      subsection[[n.subsection]] = list(title = "Multiplicity adjustment", item = item)
    }
  }

  section[[3]] = list(title = "Analysis model", subsection = subsection)

  # Section : Simulation results
  ##############################

  # Empty subsection list
  subsection = list()

  # Empty subsubsection list
  subsubsection = list()

  # Empty subsubsubsection list
  subusbsubsection = list()

  n.section = nrow(result.structure$section)
  if (!is.null(result.structure$subsection)) n.subsection = nrow(result.structure$subsection)
  else n.subsection = 0

  # Get the names of the columns to span
  span = colnames(result.structure$table.structure[[1]]$results)[which(!(colnames(result.structure$table.structure[[1]]$results) %in% c("Criterion","Test/Statistic","Result")))]

  # Create each section
  for (section.ind in 1:n.section){
    table.result.section = result.structure$table.structure[unlist(lapply(result.structure$table.structure, function(x,ind.section=section.ind) {(x$section$number == ind.section) } ))]
    # Empty subsection list
    subsection = list()
    if (n.subsection >0) {
      for (subsection.ind in 1:n.subsection){
        # Result
        item1 = list(label = "Results summary",
                     value =  table.result.section[[subsection.ind]]$results,
                     type = "table",
                     param = list(span.columns = span)
        )

        # Create a subsection
        subsection[[subsection.ind]] = list(title = table.result.section[[subsection.ind]]$subsection$title, item = list(item1))
      }
      section[[3+section.ind]] = list(title = table.result.section[[subsection.ind]]$section$title, subsection = subsection)
    }
    else {
      # Result
      item1 = list(label = "Results summary",
                   value =  table.result.section[[1]]$results,
                   type = "table",
                   param = list(span.columns = span)
      )
      subsection[[1]] = list(title = NA, item = list(item1))
      section[[3+section.ind]] = list(title = table.result.section[[1]]$section$title, subsection = subsection)
    }
  }

  # Include all sections in the report -- the report object is finalized
  report = list(title = "Clinical Scenario Evaluation", section = section)

  return(list(result.structure = result.structure, report.structure = report ))

}
# End of CreateReportStructure