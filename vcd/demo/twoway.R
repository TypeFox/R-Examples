  #####################
  ## Fourfold tables ##
  #####################

  ### Berkeley Admission Data ###
  ###############################
  data(UCBAdmissions)

  ## unstratified
  ### no margin is standardized
  x <- margin.table(UCBAdmissions, 2:1)
  fourfold(x, std = "i", extended = FALSE)
  ### std. for gender
  fourfold(x, margin = 1, extended = FALSE)
  ### std. for both
  fourfold(x, extended = FALSE)

  ## stratified
  fourfold(UCBAdmissions, extended = FALSE)
  fourfold(UCBAdmissions) ## extended plots

  ## using cotabplot
  cotabplot(UCBAdmissions, panel = function(x, condlevels, ...)
              fourfold(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
                       newpage = F, return_grob = FALSE, ...)
            )

  ### Coal Miners Lung Data ###
  #############################
  data(CoalMiners)

  ## Fourfold display, both margins equated
  fourfold(CoalMiners, mfcol = c(3,3))

  ## Log Odds Ratio Plot
  data(CoalMiners, package = "vcd")
  lor_CM <- loddsratio(CoalMiners)
  plot(lor_CM)
  lor_CM_df <- as.data.frame(lor_CM)

  # fit linear models using WLS
  age <- seq(20, 60, by = 5)
  lmod <- lm(LOR ~ age, weights = 1 / ASE^2, data = lor_CM_df)
  grid.lines(age, fitted(lmod), gp = gpar(col = "blue"))
  qmod <- lm(LOR ~ poly(age, 2), weights = 1 / ASE^2, data = lor_CM_df)
  grid.lines(age, fitted(qmod), gp = gpar(col = "red"))

  ## Fourfold display, strata equated
  fourfold(CoalMiners, std = "ind.max", mfcol = c(3,3))

  ####################
  ## Sieve Diagrams ##
  ####################

  ### Hair Eye Color ###
  ######################
  data(HairEyeColor)

  ## aggregate over `sex':
  (tab <- margin.table(HairEyeColor, 1:2))
  ## plot expected values:
  sieve(t(tab), sievetype = "expected", shade = TRUE)

  ## plot sieve diagram:
  sieve(t(tab), shade = TRUE)

  ### Visual Acuity ###
  #####################
  data(VisualAcuity)
  attach(VisualAcuity)
  sieve(Freq ~ right + left,
        data = VisualAcuity,
        subset = gender == "female",
        main = "Unaided distant vision data",
        labeling_args = list(set_varnames = c(left = "Left Eye Grade",
                               right = "Right Eye Grade")),
        shade = TRUE
        )
  detach(VisualAcuity)

  ### Berkeley Admission ###
  ##########################

  ## -> Larger tables: e.g., Cross factors
  ### Cross Gender and Admission
  data(UCBAdmissions)

  (tab <- xtabs(Freq ~ Dept + I(Gender : Admit), data = UCBAdmissions))
  sieve(tab,
        labeling_args = list(set_varnames = c("I(Gender:Admit)" = "Gender:Admission",
                               Dept = "Department")),
        main = "Berkeley Admissions Data",
        shade = TRUE
        )

  ## or use extended sieve plots:
  sieve(UCBAdmissions, shade = TRUE)

  ######################
  ## Association Plot ##
  ######################

  ### Hair Eye Color ###
  ######################
  data(HairEyeColor)
  assoc(margin.table(HairEyeColor, 1:2),
        labeling_args = list(set_varnames = c(Hair = "Hair Color", Eye = "Eye Color")),
        main = "Association Plot")

  ####################
  ## Agreement Plot ##
  ####################

  ### Sexual Fun ###
  ##################
  data(SexualFun)

  ## Kappa statistics
  Kappa(SexualFun)

  ## Agreement Chart
  agreementplot(t(SexualFun), weights = 1)
  ## Partial Agreement Chart and B-Statistics
  (agreementplot(t(SexualFun),
                      xlab = "Husband's Rating",
                      ylab = "Wife's Rating",
                      main = "Husband's and Wife's Sexual Fun")
   )

  ### MS Diagnosis data ###
  #########################
  data(MSPatients)
  ## use e.g., X11(width = 12), or expand graphics device
  agreementplot(t(MSPatients[,,1]), main = "Winnipeg Patients")
  agreementplot(t(MSPatients[,,2]), main = "New Orleans Patients")

  ##################
  ## Ternary Plot ##
  ##################

  ### sample data ###
  ###################
  (x <- rbind(c(A=10,B=10,C=80),
              c(40,30,30),
              c(20,60,20)
              )
   )
  ternaryplot(x,
              cex = 2,
              col = c("black", "blue", "red"),
              coordinates = TRUE
              )

  ### Arthritis Treatment Data ###
  ################################
  data(Arthritis)

  ## Build table by crossing Treatment and Sex
  (tab <- as.table(xtabs(~ I(Sex:Treatment) + Improved, data = Arthritis)))

  ## Mark groups
  col <- c("red", "red", "blue", "blue")
  pch <- c(1, 19, 1, 19)

  ## plot
  ternaryplot(
              tab,
              col = col,
              pch = pch,
              cex = 2,
              bg = "lightgray",
              grid_color = "white",
              labels_color = "white",
              main = "Arthritits Treatment Data"
              )
  ## legend
  grid_legend(0.8, 0.7, pch, col, rownames(tab), title = "GROUP")

  ### Baseball Hitters Data ###
  #############################
  data(Hitters)
  attach(Hitters)

  colors <- c("black","red","green","blue","red","black","blue")
  pch <- substr(levels(Positions), 1, 1)
  ternaryplot(
              Hitters[,2:4],
              pch = as.character(Positions),
              col = colors[as.numeric(Positions)],
              main = "Baseball Hitters Data"
              )
  grid_legend(0.8, 0.9, pch, colors, levels(Positions), title = "POSITION(S)")

  detach(Hitters)

  ### Lifeboats on the Titanic ###
  ################################
  data(Lifeboats)
  attach(Lifeboats)

  ternaryplot(
              Lifeboats[,4:6],
              pch = ifelse(side=="Port", 1, 19),
              col = ifelse(side=="Port", "red", "blue"),
              id  = ifelse(men/total > 0.1, as.character(boat), NA),
              dimnames_position = "edge",
              dimnames = c("Men of Crew", "Men passengers", "Women and Children"),
              main = "Lifeboats on the Titanic"
              )
  grid_legend(0.8, 0.9, c(1, 19), c("red", "blue"), c("Port", "Starboard"), title = "SIDE")

  ## Load against time for Port/Starboard boats
  plot(launch, total,
       pch = ifelse(side == "Port", 1, 19),
       col = ifelse(side == "Port", "red", "darkblue"),
       xlab = "Launch Time",
       ylab = "Total loaded",
       main = "Lifeboats on the Titanic"
       )
  legend(as.POSIXct("1912-04-15 01:48:00"), 70,
         legend = c("SIDE","Port","Starboard"),
         pch = c(NA, 1, 19),
         col = c(NA, "red", "darkblue")
         )
  text(as.POSIXct(launch),
       total,
       labels = as.character(boat),
       pos = 3,
       offset = 0.3
       )
  abline(lm(total ~ as.POSIXct(launch),
            subset = side == "Port"),
         col = "red")
  abline(lm(total ~ as.POSIXct(launch),
            subset = side == "Starboard"),
         col = "darkblue")

  detach(Lifeboats)

