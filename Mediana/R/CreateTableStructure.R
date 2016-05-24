######################################################################################################################

# Function: CreateTableStructure.
# Argument: Results returned by the CSE function and presentation model.
# Description: This function is used to create the tables for each section/subsection

CreateTableStructure = function(results = NULL, presentation.model = NULL, custom.label.sample.size = NULL, custom.label.design.parameter = NULL, custom.label.outcome.parameter = NULL, custom.label.multiplicity.adjustment = NULL){

  # TO DO: Add checks on parameters

  # Get analysis scenario grid and add a scenario number to the dataframe
  analysis.scenario.grid = results$analysis.scenario.grid
  analysis.scenario.grid$scenario = as.numeric(rownames(results$analysis.scenario.grid))
  analysis.scenario.grid$all = 1

  # Apply the label on the scenario grid
  analysis.scenario.grid.label = analysis.scenario.grid
  analysis.scenario.grid.label$design.parameter.label = as.factor(analysis.scenario.grid.label$design.parameter)
  analysis.scenario.grid.label$outcome.parameter.label = as.factor(analysis.scenario.grid.label$outcome.parameter)
  analysis.scenario.grid.label$sample.size.label = as.factor(analysis.scenario.grid.label$sample.size)
  analysis.scenario.grid.label$multiplicity.adjustment.label = as.factor(analysis.scenario.grid.label$multiplicity.adjustment)
  levels(analysis.scenario.grid.label$design.parameter.label) = custom.label.design.parameter$label
  levels(analysis.scenario.grid.label$outcome.parameter.label) = custom.label.outcome.parameter$label
  levels(analysis.scenario.grid.label$sample.size.label) = custom.label.sample.size$label
  levels(analysis.scenario.grid.label$multiplicity.adjustment.label) = custom.label.multiplicity.adjustment$label
  analysis.scenario.grid.label$all.label = as.factor(analysis.scenario.grid.label$all)

  # Create the summary table with all results
  #summary.table = CreateSummaryTable(results$evaluation.set$analysis.scenario)
  summary.table = merge(results$simulation.results,analysis.scenario.grid.label, by = c("sample.size", "outcome.parameter", "design.parameter", "multiplicity.adjustment"))
  summary.table$result = format(round(summary.table$result, digits = 4), digits = 4, nsmall = 4)

  # Check if Sample Size or event has been used to set the column names
  sample.size = (!any(is.na(results$data.structure$sample.size.set)))
  event = (!any(is.na(results$data.structure$event.set)))

  # Get the "by"
  section.by = presentation.model$section.by$by
  if (is.null(section.by)) {
    section.by = "all"
    custom.label.all = list(label = "Results", custom = FALSE)
  }
  subsection.by = presentation.model$subsection.by$by
  if (any(section.by %in% subsection.by)) stop("PresentationModel: the parameters must be defined either in the Section or in the Subsection object, but not in both")
  table.by = presentation.model$table.by$by
  if (any(section.by %in% table.by)) stop("PresentationModel: the parameters must be defined either in the Section or in the Table object, but not in both")
  if (any(subsection.by %in% table.by)) stop("PresentationModel: the parameters must be defined either in the Subsection or in the Table object, but not in both")

  # If the user used event, the "by" "event" need to be changed by "sample.size" as the
  if (event) {
    if (!is.null(section.by)) section.by = gsub("event","sample.size",section.by)
    if (!is.null(subsection.by)) subsection.by = gsub("event","sample.size",subsection.by)
    if (!is.null(table.by)) table.by = gsub("event","sample.size",table.by)
  }

  # Create a list with scenario number for all sections
  # This list will get the number of scenarios defined by the user for each parameters
  section.par = list()
  if (!is.null(section.by)){
    for (i in 1:length(section.by)){
      section.par[[i]] = levels(analysis.scenario.grid.label[,paste0(section.by[i],".label")])
    }
  } else section.par = NULL


  # Create the combination of scenario for each section
  section.create = rev(expand.grid(rev(section.par)))
  colnames(section.create) = section.by
  section.create$section = 1:nrow(section.create)
  # Create the title for the section
  section.by.label = sapply(gsub("."," ",section.by,fixed = TRUE),capwords)
  if(get(paste0("custom.label.",section.by)[1])$custom) {
    title = paste0(section.by.label[1]," (",section.create[,1],")")
  } else {
    title = paste0(section.by.label[1]," ",1:max(analysis.scenario.grid[,section.by[[1]][1]]))
  }
  if (length(section.by)>1) {
    for (i in 2:length(section.by)){
      if(get(paste0("custom.label.",section.by)[i])$custom) {
        title = paste0(title," and ",section.by.label[i]," (",section.create[,i],")")
      } else {
        title = paste0(title," and ",section.by.label[i]," ",1:max(analysis.scenario.grid[,section.by[[i]][1]]))
      }
    }
  }
  section.create$title.section = title
  if (any(section.by == "all"))  section.create$title.section = "Results"


  # Create a list with scenario number for all subsections
  subsection.create = NULL
  if (!is.null(subsection.by)){
    # This list will get the number of scenarios defined by the user for each parameters
    subsection.par = list()
    for (i in 1:length(subsection.by)){
      subsection.par[[i]] = levels(analysis.scenario.grid.label[,paste0(subsection.by[i],".label")])
    }

    # Create the combination of scenario for each subsection
    subsection.create = rev(expand.grid(rev(subsection.par)))
    colnames(subsection.create) = subsection.by
    subsection.create$subsection = 1:nrow(subsection.create)
    # Create the title for the subsection
    subsection.by.label = sapply(gsub("."," ",subsection.by,fixed = TRUE),capwords)
    if(get(paste0("custom.label.",subsection.by)[1])$custom) {
      title = paste0(subsection.by.label[1]," (",subsection.create[,1],")")
    } else {
      title = paste0(subsection.by.label[1]," ",1:max(analysis.scenario.grid[,subsection.by[[1]][1]]))
    }
    if (length(subsection.by)>1) {
      for (i in 2:length(subsection.by)){
        if(get(paste0("custom.label.",subsection.by)[i])$custom) {
          title = paste0(title," and ",subsection.by.label[i]," (",subsection.create[,i],")")
        } else {
          title = paste0(title," and ",subsection.by.label[i]," ",1:max(analysis.scenario.grid[,subsection.by[[i]][1]]))
        }
      }
    }
    subsection.create$title.subsection = title
  }

  # Create a list with order tables
  # If the user did not define any parameter to sort the table, the parameters not defined in the section or subsection will be used to sort the table by default
  table.by=c(table.by, colnames(analysis.scenario.grid.label[which(!(colnames(analysis.scenario.grid) %in% c(section.by, subsection.by, table.by, "scenario")))]))
  # If no design or no multiplicity adjustment have been defined, we can delete them from the table.by
  if (is.null(results$analysis.structure$mult.adjust)) table.by=table.by[which(table.by!="multiplicity.adjustment")]
  if (is.null(results$data.structure$design.parameter.set)) table.by=table.by[which(table.by!="design.parameter")]
  if (any(section.by != "all")) table.by=table.by[which(table.by!="all")]

  table.create = NULL
  if (length(table.by)>0){
    # This list will get the number of scenarios defined by the user for each parameters
    table.par = list()
    for (i in 1:length(table.by)){
      table.par[[i]] = as.numeric(levels(as.factor(analysis.scenario.grid.label[,table.by[i]])))
    }
    # Create the combination of scenario for each table
    table.create = rev(expand.grid(rev(table.par)))
    colnames(table.create) = table.by
  }


  # Create report structure
  if(!is.null(subsection.create)){
    report.structure = rev(merge(subsection.create,section.create,by=NULL,suffixes = c(".subsection",".section")))
  } else report.structure = section.create

  # Get the scenario number for each section/subsection
  report.structure.scenario = list()
  for (i in 1:nrow(report.structure)){
    report.structure.temp = as.data.frame(report.structure[i,])
    colnames(report.structure.temp) = paste0(colnames(report.structure),".label")
    report.structure.scenario[[i]] = as.vector(merge(analysis.scenario.grid.label,report.structure.temp)$scenario)
  }

  # Create a list containing the table to report under each section/subsection
  report.structure.scenario.summary.table = list()
  for (i in 1:nrow(report.structure)){
    report.structure.scenario.summary.table[[i]] = summary.table[which(summary.table$scenario %in% report.structure.scenario[[i]]),c(table.by,"criterion","test.statistic","result")]
    colnames(report.structure.scenario.summary.table[[i]])=c(sapply(gsub("."," ",table.by,fixed = TRUE),capwords),"Criterion","Test/Statistic","Result")
    rownames(report.structure.scenario.summary.table[[i]])=NULL
  }

  # Order the table as requested by the user
  if (length(table.by)>0){
    table.by.label = sapply(gsub("."," ",table.by,fixed = TRUE),capwords)
    data.order = as.data.frame(report.structure.scenario.summary.table[[1]][,table.by.label])
    colnames(data.order) = table.by
    order.table = order(apply(data.order, 1, paste, collapse = ""))
    report.structure.scenario.summary.table.order = lapply(report.structure.scenario.summary.table,function(x) x[order.table,])
  } else report.structure.scenario.summary.table.order = report.structure.scenario.summary.table

  # Add the labels
  report.structure.scenario.summary.table.order = lapply(report.structure.scenario.summary.table.order, function(x) {
    if ("Design Parameter" %in% colnames(x)) {
      x[,"Design Parameter"] = as.factor(x[,"Design Parameter"])
      levels(x[,"Design Parameter"]) = custom.label.design.parameter$label
    }
    if ("Outcome Parameter" %in% colnames(x)) {
      x[,"Outcome Parameter"] = as.factor(x[,"Outcome Parameter"])
      levels(x[,"Outcome Parameter"]) = custom.label.outcome.parameter$label
    }
    if ("Sample Size" %in% colnames(x)) {
      x[,"Sample Size"] = as.factor(x[,"Sample Size"])
      levels(x[,"Sample Size"]) = custom.label.sample.size$label
    }
    if ("Multiplicity Adjustment" %in% colnames(x)) {
      x[,"Multiplicity Adjustment"] = as.factor(x[,"Multiplicity Adjustment"])
      levels(x[,"Multiplicity Adjustment"]) = custom.label.multiplicity.adjustment$label
    }

    return(x) })


  # Change the Sample Size column name if Event has been used
  ChangeColNames = function(x) {
    colnames(x)[colnames(x)=="Sample Size"] <- "Event Set"
    x
  }
  if (event) report.structure.scenario.summary.table.order = lapply(report.structure.scenario.summary.table.order, ChangeColNames)

  # Create the object to return, i.e. a list with the parameter of the section/subsection and the table
  res = list()
  for (i in 1:nrow(report.structure)){
    report.structure.temp = as.data.frame(report.structure[i,])
    colnames(report.structure.temp) = colnames(report.structure)
    res[[i]] = list(section = list(number = report.structure.temp$section, title = report.structure.temp$title.section),
                    subsection = list(number = report.structure.temp$subsection, title = report.structure.temp$title.subsection),
                    parameter = as.character(report.structure.temp[,c(section.by, subsection.by)]),
                    results = report.structure.scenario.summary.table.order[[i]])
  }

  # Return the list of results
  return(list(section = section.create, subsection = subsection.create, table.structure = res))
}
# End of CreateTableStructure
