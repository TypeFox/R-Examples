dataprep <-
function(
                        foo = NULL,
                        predictors = NULL,
                        predictors.op = "mean",
                        special.predictors = NULL,
                        dependent = NULL,
                        unit.variable = NULL,
                        time.variable = NULL,
                        treatment.identifier = NULL,
                        controls.identifier = NULL,
                        time.predictors.prior = NULL,
                        time.optimize.ssr = NULL,
                        time.plot = time.optimize.ssr,
                        unit.names.variable = NA
                        )
  
  {

    if(is.data.frame(foo) == FALSE ){stop("\n No data.frame supplied in foo.\n")}

    # identify unit.variable
    if(mode(unit.variable) == "character"){unit.variable<-which(names(foo)==unit.variable)}
    if(is.null(unit.variable) == TRUE ||  mode(foo[,unit.variable])!="numeric")
      { stop("\n unit.variable not found as numeric variable in foo.\n")}
    if(length(unit.variable)!=1){stop("\ Please specify only one unit.variable\n")}

    # identify time.variable
    if(mode(time.variable) == "character"){time.variable<-which(names(foo)==time.variable)}
    if(is.null(time.variable) == TRUE || mode(foo[,time.variable])!="numeric" )
      {
       stop("\n time.variable not found as numeric variable in foo.\n")
      }
     if(length(time.variable)!=1){stop("\ Please specify only one time.variable\n")}
    # identify units
     # check if unit.name var is there (required if any identifier is given as character)
     if(
        mode(treatment.identifier) == "character" ||
        mode(controls.identifier)  == "character" ||
        is.na(unit.names.variable) == FALSE
        )
      {
       if(mode(unit.names.variable) == "character")
         {unit.names.variable<-which(names(foo)==unit.names.variable)}
       if(is.na(unit.names.variable) == TRUE || mode(foo[,unit.names.variable])!="character")
         {stop("\n unit.names.variable not found as character variable in foo.\n")}
      }
     # check treated
     if(length(treatment.identifier) != 1){stop("\n please specify a single treated unit\n")}
     if(mode(treatment.identifier) == "character")
      {
       if(treatment.identifier %in% foo[,unit.names.variable] == FALSE  )
        {stop("\n treated unit not found in unit.names.variable.\n")}
       tr.id = unique(foo[foo[,unit.names.variable] == treatment.identifier,unit.variable])
       if(length(tr.id)!=1)
         {stop(paste("\n treatment unit name",treatment.identifier,"refers to multiple numbers in unit.variable"))}
       treatment.identifier.name = treatment.identifier
       treatment.identifier = tr.id
      } else {
       if(treatment.identifier %in% foo[,unit.variable] == FALSE )
        {stop("\n treated unit not found in unit.variable.\n")}
        if(is.na(unit.names.variable) == TRUE )
         {
           treatment.identifier.name = treatment.identifier
          } else {
           treatment.identifier.name = unique(foo[treatment.identifier==foo[,unit.variable],unit.names.variable])
           if(length(treatment.identifier.name)>1)
            {stop("\n treatment.identifier refers to multiple names in unit.names.variable")}
          }
      }

     # check controls
     if(length(controls.identifier) < 2 ){stop("\n please specify at least two control units\n")}
     if(sum(duplicated(controls.identifier))>0)
      {stop("\n duplicate control units in controls.identifier\n")}
     if(mode(controls.identifier) == "character")
       { co.store <- c()
         for(i in controls.identifier)
          {
            if(i %in% foo[,unit.names.variable] == FALSE )
             {stop(paste("\n control unit",i,"not found in unit.names.variable"))}
            co.id = unique(foo[foo[,unit.names.variable] == i,unit.variable])
            if(length(co.id) != 1)
             {stop(paste("\n control unit name",i," refers to multiple numbers in unit.variable"))}
            co.store <- c(co.store,co.id)
          }
         controls.identifier.name = controls.identifier
         controls.identifier = co.store
        } else {
         co.store <- c()
         for(i in controls.identifier)
          {
            if(i %in% foo[,unit.variable] == FALSE )
             {stop(paste("\n control unit",i,"not found in unit.variable"))}
            if(is.na(unit.names.variable) == FALSE )
            {
             co.id = unique(foo[foo[,unit.variable] == i,unit.names.variable])
             if(length(co.id) != 1)
              {stop(paste("\n control unit number",i," refers to multiple names in unit.names.variable"))}
             co.store <- c(co.store,co.id)
             }
          }
         if(is.na(unit.names.variable) == FALSE )
          {
           controls.identifier.name = co.store
         } else {
           controls.identifier.name = controls.identifier
         }
        }

        # more checks
        if(treatment.identifier.name %in% controls.identifier.name)
         {stop("\n treated unit among controls\n")}
         
        if(sum(duplicated(c(controls.identifier.name,treatment.identifier.name))) > 0)
         {stop("n duplicate unit.variable.names across units\n")}
         
        # sort first by unit, then by time variable
         foo[,time.variable] <- as.numeric(as.character(foo[,time.variable]))
         foo[,unit.variable] <- as.numeric(as.character(foo[,unit.variable]))
         foo <- foo[order(foo[,unit.variable], foo[,time.variable]),]

      # Get rows
       treatment.rows <- which(foo[,unit.variable] %in% treatment.identifier)
       control.rows   <- which(foo[,unit.variable] %in% controls.identifier)

      # check if panel is unbalanced  (possible cases where this too restrictive, but imposes discipline)
       balcheck <-       table(    foo[c(control.rows,treatment.rows),unit.variable],
                          foo[c(control.rows,treatment.rows),time.variable])         
        if( length(unique(balcheck)) != 1 || unique(balcheck)!= 1)
      {
        stop("\n Your panel, as described by unit.variable and time.variable, is unbalanced. Balance it and run again.\n")
      }

     # Now check and get time identifiers
      if(sum(is.null(time.predictors.prior))> 0)
       {stop("time.predictors.prior missing")}
      if(sum(is.null(time.optimize.ssr))>0)
       {stop("time.optimize.ssr missing")}
      if(sum(is.null(time.plot))>0)
       {stop("time.plot missing")}
       
      t.list <- list(time.predictors.prior,time.optimize.ssr,time.plot)
      names(t.list) <- c("predictors.prior","optimize.ssr","plot")
      for(i in 1:3)
       {
        if(mode(t.list[[i]]) != "numeric")
         {stop(paste("\n time",names(t.list[i]),"must be numeric\n"))}
        if(sum(duplicated(t.list[[i]]))>0)
         {stop(paste("\n duplicates in time",names(t.list[i]),"\n"))}
        if(length(t.list[[i]]) < 1)
         {stop(paste("\n specificy at least one period in time",names(t.list[i]),"\n"))}
        for(p in t.list[[i]])
         {
          if(p %in% unique(foo[,time.variable]) == FALSE)
           stop(paste("\n time period ",p," from time.",names(t.list[i])," not found in time.variable\n",sep=""))
         }
        }
     
     # get time rows
     time.predictors.prior.rows <-
      which(foo[,time.variable] %in% time.predictors.prior)

     time.optimize.ssr.rows <-
      which(foo[,time.variable] %in% time.optimize.ssr)

    # Build X1 & X0 Predictor matrices
    
   # first, run check on regular predictors, if specified
   if(is.null(predictors)==FALSE)
   {
     # predictor checks
     if(sum(duplicated(predictors))>0)
      {stop("\n duplicates found in predictors \n")}

     for(i in predictors)
      {
        if(mode(foo[,i]) != "numeric")
         {stop(paste("\n predictor",i,"not found as numeric variable in foo \n"))}
      }

      if(mode(predictors)=="character")
       {
         pred.no <- c()
         for(i in 1:length(predictors)){pred.no <- c(pred.no,which(names(foo) == predictors[i]))}
         predictors <- pred.no
       }
    # X1 matrix for treated
    X1 <-
      data.frame(foo[
                        intersect(
                                  treatment.rows,
                                  time.predictors.prior.rows
                                  ),
                        predictors
                        ]
                    )
    colnames(X1) <- names(foo)[predictors]
    rownames(X1) <- time.predictors.prior

    # deep missing checker:
    for(i in 1:ncol(X1))
     {
       if(sum(is.na(X1[,i])) == nrow(X1))
        {stop(paste("\n predictor",names(X1)[i],"has missing data for all periods in time.predictors.prior\n"))}
        
      for(j in 1:nrow(X1))
       {
         if(is.na(X1[j,i])){
         cat(paste("\n Missing data- treated unit; predictor:",names(X1)[i],"; for period:",rownames(X1)[j],
         "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
       }
     }

    # aggregate
    X1 <- apply(X1, 2, paste(predictors.op), na.rm = TRUE)
    X1 <- as.matrix(X1)
    rownames(X1) <- names(foo)[predictors]
    colnames(X1) <- treatment.identifier
    
    } else { # if no regular predictors are specified, pass a void matrix to special predictors
    
    X1 <- matrix(NA,0,1)
    colnames(X1) <- treatment.identifier 
    }
    
    # X0 matrix for controls
    if(is.null(predictors)==FALSE)
    {
    X0 <- data.frame(foo[intersect(control.rows,
                                      time.predictors.prior.rows
                                      ),
                            c(predictors, unit.variable)
                            ]
                        )
     checknames <- rep(time.predictors.prior,length(controls.identifier))

     # deep missing checker:
    for(i in 1:(ncol(X0)-1))
     {

       for(p in controls.identifier)
        {
         if(sum(is.na(X0[X0[,ncol(X0)] == p,i])) == length(is.na(X0[X0[,ncol(X0)] == p,i])))
         {stop(paste("\n control unit:",p,"; predictor:",names(X0)[i],"has missing data for all periods in time.predictors.prior\n"))}
        }
       for(j in 1:nrow(X1))
       {
         if(is.na(X0[j,i])){
         cat(paste("\n Missing data - control unit:",X0[j,ncol(X0)],"; predictor:",names(X0)[i],"; for period:",checknames[j],
         "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
       }
     }

    X0 <- split(X0, X0[,dim(X0)[2]])
    X0 <- sapply(X0, apply, 2, mean, na.rm = TRUE, simplify = TRUE)
    X0 <- as.matrix(X0[-dim(X0)[1],])
    
    
   # Take transpose in presence of a single predictor only
   if(dim(X0)[2]==1)
   {
    X0 <- t(X0)
    rownames(X0) <- names(foo)[predictors]
    colnames(X0) <- controls.identifier
   }

    } else { # if no regular predictors are specified, pass a void matrix to special predictors
    X0 <- matrix(NA,0,length(controls.identifier))
    #rownames(X0) <- names(foo)[predictors]
    colnames(X0) <- controls.identifier
    }

    # Add Special Predictors to X1 and X0

    if(is.null(special.predictors) == FALSE)
      {
        if(is.list(special.predictors) == FALSE)
         {stop("\nspecial.predictors is not a list object\n")}

        for(i in 1:length(special.predictors))
          {
          
           # checks for special predictors
            if(is.list(special.predictors[[i]]) == FALSE)
             {stop(paste("\n special.predictor number",i,"is not a list object\n"))}

            if(length(special.predictors[[i]]) != 3)
             {stop(paste("\n special.predictor number",i,"is misspecified (requires predictor name, time periods, and name of operation\n"))}

           # name check
            sp.name <- special.predictors[[i]][[1]]
            if(is.na(sp.name) || length(sp.name) != 1)
            {stop(paste("\n predictor name",sp.name,"of special.predictor number",i,"misspecified\n"))}
            
           # availability check
           if(mode(foo[,sp.name]) != "numeric")
            {stop(paste("\n special predictor named ",sp.name,"not found as numeric variable in foo \n"))}
           
           # time check
            sp.time <- special.predictors[[i]][[2]]
            if(mode(sp.time) != "numeric")
            {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") the time period is misspecified\n"))}
            if(sum(duplicated(sp.time))>0)
             {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") time period contains duplicates\n"))}
            if(length(sp.time)<1)
             {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") specify at least on time period\n"))}
             
           # time availability check
           for(p in sp.time)
            {
             if(p %in% unique(foo[,time.variable]) == FALSE)
              {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") time period",p,"not found in time.variable\n"))}
            }

           # operator check
           sp.op <- special.predictors[[i]][[3]]
            if(mode(sp.op) != "character" || length(sp.op) !=1 )
            {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") the operator is mispecified\n"))}

           # now go and built matrices
           spf <- spec.pred.func(list.object = special.predictors[[i]],
                                  tr.numb = treatment.identifier,
                                  co.numb = controls.identifier,
                                  unit.var = unit.variable,
                                  time.var = time.variable,
                                  foo.object = foo,
                                  X0.inner = X0,
                                  X1.inner = X1
                                  )

            X0 <- spf[[1]]
            X1 <- spf[[2]]

          }

      }

   # no predictors check
    if(nrow(X0)==0)
     {stop("No predictors specified. Please specify at least on predictor")}

   # Dependent Variable Pretreatment

    if(is.null(dependent))
     {stop("\n dependent variabale is missing")}

    if(mode(foo[,dependent]) != "numeric")
         {stop(paste("\n dependent variable",dependent,"not found as numeric variable in foo \n"))}

    if(mode(dependent) == "character")
      {dependent <- which(is.element(names(foo), dependent))}
    #foo[,dependent] <- as.numeric(as.character(foo[,dependent]))

    Z1 <- foo[
              intersect(
                        treatment.rows,
                        time.optimize.ssr.rows
                        ),
              dependent
              ]

    Z1 <- as.matrix(Z1)

    # missing check
    if(sum(is.na(Z1)) == length(Z1))
     {stop("\n treated unit: dependent variable is missing for all periods in time.optimize.ssr.rows \n")}

    rownames(Z1) <- time.optimize.ssr
    colnames(Z1) <- treatment.identifier
    
    for(i in 1:length(Z1))
      {
        if(is.na(Z1[i]))
        {stop(paste("\n treated unit: dependent variable is missing for period",rownames(Z1)[i],"in time.optimize.ssr.rows \n"))}
      }

    Z0 <-
      as.matrix(foo[
                    intersect(
                              control.rows,
                              time.optimize.ssr.rows
                              ),
                    c(dependent)
                    ]
                )

    Z0 <- matrix(
                 Z0,
                 byrow = FALSE,
                 nrow = length(time.optimize.ssr),
                 ncol = length(controls.identifier),
                 )
                 
    colnames(Z0)<- controls.identifier
    rownames(Z0)<- time.optimize.ssr

    for(i in 1:ncol(Z0))
      {
        for(j in 1:nrow(Z0))
        {
        if(is.na(Z0[j,i]))
        {stop(paste("\n control unit:",colnames(Z0)[i],"; dependent variable",dependent,"is missing for period",rownames(Z0)[j],"in time.optimize.ssr.rows \n"))}
        }
      }
    

     # dependent Variable Plot
        time.plot.rows <-
          which(is.element(foo[,time.variable], time.plot ) == TRUE)

        Y1plot <- foo[intersect(treatment.rows, time.plot.rows), dependent]
        Y1plot <- as.matrix(Y1plot)
        rownames(Y1plot) <- time.plot
        colnames(Y1plot) <- treatment.identifier

     if(sum(is.na(Y1plot)) == length(Y1plot))
     {stop("\n treated unit: dependent variable is missing for all periods in time.plot \n")}

      for(i in 1:length(Y1plot))
       {
        if(is.na(Y1plot[i]))
        {stop(paste("\n treated unit: dependent variable is missing for period",rownames(Y1plot)[i],"in time.plot \n"))}
       }

        Y0plot <- as.matrix(foo[
                                intersect(control.rows, time.plot.rows),
                                c(dependent)
                                ]
                            )

        Y0plot  <- matrix(
                          Y0plot,
                          byrow = FALSE,
                          nrow = length(time.plot),
                          ncol = length(controls.identifier),
                          )


        rownames(Y0plot) <- time.plot
        colnames(Y0plot) <- controls.identifier
     

        # table with unit names
              names.and.numbers <-
        data.frame(c(treatment.identifier.name,controls.identifier.name),
                   c(treatment.identifier,controls.identifier))



      names(names.and.numbers) <- c(
                                    "unit.names",
                                    "unit.numbers"
                                    )



  #######################
       
  tag <- list(
              foo = as.character(foo),
              predictors = predictors,
              predictors.op = predictors.op,
              special.predictors = special.predictors,
              dependent = dependent,
              unit.variable = unit.variable,
              time.variable = time.variable,
              treatment.identifier = treatment.identifier,
              controls.identifier = controls.identifier,
              time.predictors.prior = time.predictors.prior,
              time.optimize.ssr = time.optimize.ssr,
              time.plot = time.plot,
              unit.names.variable = unit.names.variable
              )

######################################

  output <- list(
                 X0 = X0,
                 X1 = X1,
                 Z0 = Z0,
                 Z1 = Z1,
                 Y0plot = Y0plot,
                 Y1plot = Y1plot,
                 names.and.numbers = names.and.numbers,
                 tag = tag
                 )

  return(invisible(output))

}

