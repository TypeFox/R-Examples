if (getRversion() >= "2.15.1") utils::globalVariables(c("counter", "Prediction",
                                                        "input.data", "old.d"))

DynNom.glm <- function(model, data,
                       clevel = 0.95, covariate = c("slider", "numeric")) {

  data <- data.frame(data)

  if (length(dim(data)) > 2)
    stop("Error in data format: dataframe format required")

  if (attr(model$terms, "dataClasses")[[1]] == "logical")
    stop("Error in model syntax: logical form for response not supported")

  if (tail(names(attr(model$terms,"dataClasses")),n=1)=="(weights)") {
    n.terms <- length(attr(model$terms,"dataClasses"))
    attr(model$terms,"dataClasses") <- attr(model$terms,"dataClasses")[1:n.terms - 1]
  }

  for(i in 1:length(names(attr(model$terms, "dataClasses")))) {
    com1=numeric(length(names(data)))
    for(j in 1:length(names(data))) {
      if (names(attr(model$terms, "dataClasses"))[i]==names(data)[j]) com1[j]=1
    }
    if (sum(com1)==0)
      stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
  }

  covariate <- match.arg(covariate)
  input.data <- NULL
  old.d <- NULL

  runApp(list(

    ui = bootstrapPage(fluidPage(
      titlePanel("Dynamic Nomogram"),
      sidebarLayout(sidebarPanel(uiOutput("manySliders.f"),
                                 uiOutput("manySliders.n"),
                                 checkboxInput("limits", "Set x-axis ranges"),
                                 conditionalPanel(condition = "input.limits == true",
                                                  numericInput("uxlim", "x-axis lower", NA),
                                                  numericInput("lxlim", "x-axis upper", NA)),
                                 actionButton("add", "Predict"),
                                 br(), br(),
                                 helpText("Press Quit to exit the application"),
                                 actionButton("quit", "Quit")
      ),
      mainPanel(tabsetPanel(id = "tabs",
                            tabPanel("Graphical Summary", plotOutput("plot")),
                            tabPanel("Numerical Summary", verbatimTextOutput("data.pred")),
                            tabPanel("Model Summary", verbatimTextOutput("summary"))
      )
      )
      ))),

    server = function(input, output){
      q <- observe({
        if (input$quit == 1)
          stopApp()
      })

      limits0 <- c(mean(as.numeric(model$model[[1]])) - 3 * sd(model$model[[1]]),
                   mean(as.numeric(model$model[[1]])) + 3 * sd(model$model[[1]]))
      limits <- reactive({
        if (as.numeric(input$limits) == 1) {
          limits <- c(input$lxlim, input$uxlim)
        } else {
          limits <- limits0
        }
      })

      neededVar <- names(attr(model$terms, "dataClasses"))
      data <- data[, neededVar]
      input.data <<- data[1, ]
      input.data[1, ] <<- NA

      b <- 1
      i.factor <- NULL
      i.numeric <- NULL
      for (j in 2:length(attr(model$terms, "dataClasses"))) {
        for (i in 1:length(data)) {
          if (names(attr(model$terms, "dataClasses"))[j] == names(data)[i]) {
            if (attr(model$terms, "dataClasses")[[j]] == "factor" |
                attr(model$terms, "dataClasses")[[j]] == "ordered" |
                attr(model$terms, "dataClasses")[[j]] == "logical") {
              i.factor <- rbind(i.factor, c(names(attr(model$terms, "dataClasses"))[j], j, i, b))
              (break)()
            }
            if (attr(model$terms, "dataClasses")[[j]] == "numeric") {
              i.numeric <- rbind(i.numeric, c(names(attr(model$terms, "dataClasses"))[j], j, i))
              b <- b + 1
              (break)()
            }
          }
        }
      }

      nn <- nrow(i.numeric)
      if (is.null(nn)) {
        nn <- 0
      }
      nf <- nrow(i.factor)
      if (is.null(nf)) {
        nf <- 0
      }

      if (nf > 0) {
        output$manySliders.f <- renderUI({
          slide.bars <- list(lapply(1:nf, function(j) {
            selectInput(paste("factor", j, sep = ""),
                        names(attr(model$terms, "dataClasses")[as.numeric(i.factor[j, 2])]),
                        model$xlevels[[as.numeric(i.factor[j, 2]) - as.numeric(i.factor[j, 4])]], multiple = FALSE)
          }))
          do.call(tagList, slide.bars)
        })
      }

      if (nn > 0) {
        output$manySliders.n <- renderUI({
          if (covariate == "slider") {
            slide.bars <- list(lapply(1:nn, function(j) {
              sliderInput(paste("numeric", j, sep = ""),
                          names(attr(model$terms, "dataClasses")[as.numeric(i.numeric[j, 2])]),
                          min = as.integer(min(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
                          max = as.integer(max(na.omit(data[, as.numeric(i.numeric[j, 3])]))) + 1,
                          value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
            }))
          }
          if (covariate == "numeric") {
            slide.bars <- list(lapply(1:nn, function(j) {
              numericInput(paste("numeric", j, sep = ""),
                           names(attr(model$terms, "dataClasses")[as.numeric(i.numeric[j, 2])]),
                           value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
            }))
          }
          do.call(tagList, slide.bars)
        })
      }

      a <- 0
      new.d <- reactive({
        if (nf > 0) {
          input.f <- vector("list", nf)
          for (i in 1:nf) {
            input.f[[i]] <- local({
              input[[paste("factor", i, sep = "")]]
            })
            names(input.f)[i] <- i.factor[i, 1]
          }
        }
        if (nn > 0) {
          input.n <- vector("list", nn)
          for (i in 1:nn) {
            input.n[[i]] <- local({
              input[[paste("numeric", i, sep = "")]]
            })
            names(input.n)[i] <- i.numeric[i, 1]
          }
        }
        if (nn == 0) {
          out <- data.frame(do.call("cbind", input.f))
        }
        if (nf == 0) {
          out <- data.frame(do.call("cbind", input.n))
        }
        if (nf > 0 & nn > 0) {
          out <- data.frame(do.call("cbind", input.f), do.call("cbind", input.n))
        }
        if (a == 0) {
          wher <- match(names(out), names(input.data)[-1])
          out <- out[wher]
          input.data <<- rbind(input.data[-1], out)
        }
        if (a > 0) {
          wher <- match(names(out), names(input.data))
          out <- out[wher]
          input.data <<- rbind(input.data, out)
        }
        a <<- a + 1
        out
      })

      p1 <- NULL
      old.d <- NULL
      data2 <- reactive({
        if (input$add == 0)
          return(NULL)
        if (input$add > 0) {
          if (isTRUE(compare(old.d, new.d())) == FALSE) {
            OUT <- isolate({
              pred <- predict(model, newdata = new.d(), type = "link", conf.int = clevel, se.fit = TRUE)
              lwb <- pred$fit - (qnorm(1 - (1 - clevel)/2) * pred$se.fit)
              upb <- pred$fit + (qnorm(1 - (1 - clevel)/2) * pred$se.fit)
              d.p <- data.frame(Prediction = model$family$linkinv(pred$fit),
                                Lower.bound = model$family$linkinv(lwb), Upper.bound = model$family$linkinv(upb))
              old.d <<- new.d()
              data.p <- cbind(d.p, counter = 1)
              p1 <<- rbind(p1, data.p)
              p1$count <- seq(1, dim(p1)[1])
              p1
            })
          } else {
            p1$count <- seq(1, dim(p1)[1])
            OUT <- p1
          }
        }
        OUT
      })

      if (model$family$family == "gaussian" |
          model$family$family == "inverse.gaussian" |
          model$family$family == "quasi") {
        output$plot <- renderPlot({
          if (input$add == 0)
            return(NULL)
          OUT <- isolate({
            if (is.null(new.d()))
              return(NULL)
            if (is.na(input$lxlim) | is.na(input$uxlim)) {
              lim <- limits0
            } else {
              lim <- limits()
            }
            yli <- c(0 - 0.5, 10 + 0.5)
            if (dim(input.data)[1] > 11)
              yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
            p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
            p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
            p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
            p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
                                    size = 1.45, height = 0.4, colour = data2()$count)
            p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for Response", sep = ""),
                          x = "Response", y = NULL)
            p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 14))
            print(p)
          })
        })
      }

      if (model$family$family == "poisson" |
          model$family$family == "quasipoisson" |
          model$family$family == "Gamma") {
        output$plot <- renderPlot({
          if (input$add == 0)
            return(NULL)
          OUT <- isolate({
            if (is.null(new.d()))
              return(NULL)
            if (is.na(input$lxlim) | is.na(input$uxlim)) {
              lim <- c(0, limits0[2])
            } else {
              lim <- limits()
            }
            yli <- c(0 - 0.5, 10 + 0.5)
            if (dim(input.data)[1] > 11)
              yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
            p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
            p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
            p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
            p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
                                    size = 1.45, height = 0.4, colour = data2()$count)
            p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for Response", sep = ""),
                          x = "Response", y = NULL)
            p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 14))
            print(p)
          })
        })
      }

      if (model$family$family == "binomial" |
          model$family$family == "quasibinomial") {
        output$plot <- renderPlot({
          if (input$add == 0)
            return(NULL)
          OUT <- isolate({
            if (is.null(new.d()))
              return(NULL)
            if (is.na(input$lxlim) | is.na(input$uxlim)) {
              lim <- c(0, 1)
            } else {
              lim <- limits()
            }
            yli <- c(0 - 0.5, 10 + 0.5)
            if (dim(input.data)[1] > 11)
              yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
            p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
            p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
            p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
            p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
                                    size = 1.45, height = 0.4, colour = data2()$count)
            p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for Response", sep = ""),
                          x = "Probability", y = NULL)
            p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 14))
            print(p)
          })
        })
      }

      output$data.pred <- renderPrint({
        if (input$add > 0) {
          OUT <- isolate({
            if (nrow(data2() > 0)) {
              if (dim(input.data)[2] == 1) {
                in.d <- data.frame(input.data[-1, ])
                names(in.d) <- names(attr(model$terms, "dataClasses"))[2]
                data.p <- cbind(in.d, data2()[1:3])
              }
              if (dim(input.data)[2] > 1) {
                data.p <- cbind(input.data[-1, ], data2()[1:3])
              }
              stargazer(data.p, summary = FALSE, type = "text")
            }
          })
        }
      })

      output$summary <- renderPrint({
        if (model$family$family == "binomial" |
            model$family$family == "quasibinomial") {
          coef.c <- exp(model$coef)
          ci.c <- exp(suppressMessages(confint(model, level = clevel)))
          stargazer(model, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F, type = "text",
                    omit.stat = c("LL", "ser", "f"), ci = TRUE, ci.level = clevel, single.row = TRUE,
                    title = paste(model$family$family, " regression (", model$family$link, "): ", model$formula[2],
                                  " ", model$formula[1], " ", model$formula[3], sep = ""))
        } else {
          stargazer(model, type = "text", omit.stat = c("LL", "ser", "f"),
                    ci = TRUE, ci.level = clevel, single.row = TRUE,
                    title = paste(model$family$family, " regression (", model$family$link, "): ", model$formula[2],
                                  " ", model$formula[1], " ", model$formula[3], sep = ""))
        }
      })}
  )
  )
}
