utils::globalVariables(c("counter", "Prediction", "input.data", "old.d", "old.d2", "surv", "event", "n.risk", "part"))

DynNom.coxph <- function(model, data,
                         clevel = 0.95, covariate = c("slider", "numeric"),
                         ptype = c("st", "1-st")) {

  data <- data.frame(data)

  if (length(dim(data)) > 2)
    stop("Error in data format: dataframe format required")

  if (attr(model$terms, "dataClasses")[[1]] == "logical")
    stop("Error in model syntax: logical form for response not supported")

  if (tail(names(attr(model$terms,"dataClasses")),n=1)=="(weights)") {
    n.terms <- length(attr(model$terms,"dataClasses"))
    attr(model$terms,"dataClasses") <- attr(model$terms,"dataClasses")[1:n.terms - 1]
  }

  if (attr(model, "class")[1] == "coxph.null") {
    stop("Error in model syntax: the model is null")
  }

  n.strata <- length(attr(model$terms, "specials")$strata)
  dim.terms <- length(names(attr(model$terms, "dataClasses")))

  for (i in 2:dim.terms) {
    if (substr(names(attr(model$terms, "dataClasses"))[i], 1, 6) == "strata") {
      nch <- nchar(names(attr(model$terms, "dataClasses"))[i])
      names(attr(model$terms, "dataClasses"))[i] <- substr(names(attr(model$terms,
                                                                      "dataClasses"))[i], 8, (nch - 1))
    }
  }

  if (!is.null(attr(model$terms, "specials")$tt)) {
    stop("Error in model syntax: coxph models with a time dependent covariate is not supported")
  }

  for(i in 2:length(names(attr(model$terms, "dataClasses")))) {
    com1=numeric(length(names(data)))
    for(j in 1:length(names(data))) {
      if (names(attr(model$terms, "dataClasses"))[i]==names(data)[j]) com1[j]=1
    }
    if (sum(com1)==0)
      stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
  }

  covariate <- match.arg(covariate)
  ptype <- match.arg(ptype)
  input.data <- NULL
  old.d <- NULL

  n.strata <- length(attr(model$terms, "specials")$strata)

  runApp(list(

    ui = bootstrapPage(fluidPage(
      titlePanel("Dynamic Nomogram"),
      sidebarLayout(sidebarPanel(uiOutput("manySliders.f"),
                                 uiOutput("manySliders.n"),
                                 checkboxInput("trans", "Alpha blending (transparency)", value = TRUE),
                                 actionButton("add", "Predict"),
                                 br(), br(),
                                 helpText("Press Quit to exit the application"),
                                 actionButton("quit", "Quit")
      ),
      mainPanel(tabsetPanel(id = "tabs",
                            tabPanel("Estimated S(t)", plotOutput("plot")),
                            tabPanel("Predicted Survival", plotOutput("plot2")),
                            tabPanel("Numerical Summary", verbatimTextOutput("data.pred")),
                            tabPanel("Model Summary", verbatimTextOutput("summary"))
      )
      )
      ))),

    server = function(input, output){

      observe({
        if (input$quit == 1)
          stopApp()
      })

      neededVar <- names(attr(model$terms, "dataClasses"))[-1]
      if (length(attr(model$terms, "term.labels")) == 1) {
        input.data <<- data.frame(data[1, neededVar])
        names(input.data)[1] <<- names(attr(model$terms, "dataClasses"))[-1]
      } else {
        input.data <<- data[1, neededVar]
      }
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

      tt <- names(attr(model$terms, "dataClasses"))[1]
      dd <- unlist(strsplit(substr(tt, 6, nchar(tt) - 1), "[,]"))
      tim <- dd[1]
      sts <- substr(dd[2], 2, nchar(dd[2]))

      if (length(attr(model$terms, "term.labels")) == 1) {
        input.data <<- data.frame(cbind(stt = NA, ti = NA, cov = NA), NO=NA)
        names(input.data)[3] <<- paste(attr(model$terms, "term.labels"))
        names(input.data)[1:2] <<- c(paste(sts), paste(tim))
      } else {
        data1 <- data[, neededVar]
        input.data <<- cbind(stt = NA, ti = NA, data1[1, ], NO=NA)
        names(input.data)[1:2] <<- c(paste(sts), paste(tim))
        input.data[1, ] <<- NA
      }

      if (length(i.numeric) == 0) {
        i.numeric <- matrix(ncol = 3)
        i.numeric <- rbind(i.numeric, V1 = paste(tim))
        i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1],1])
        i.numeric <- rbind(i.numeric, V1 = paste(sts))
        i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1], 1])
        i.numeric <- i.numeric[-1, ]
      } else {
        i.numeric <- rbind(i.numeric, V1 = paste(tim))
        i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1], 1])
        i.numeric <- rbind(i.numeric, V1 = paste(sts))
        i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1], 1])
      }

      limits0 <- c(0, as.integer(quantile(na.omit(data[,as.numeric(i.numeric[dim(i.numeric)[1] - 1, 3])]), probs = 0.7)))

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

      if (nn > 1) {
        output$manySliders.n <- renderUI({
          if (covariate == "slider") {
            if (nn > 2){
              slide.bars <- list(lapply(1:(nn - 2), function(j) {
                sliderInput(paste("numeric", j, sep = ""), i.numeric[j, 1],
                            min = as.integer(min(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
                            max = as.integer(max(na.omit(data[, as.numeric(i.numeric[j, 3])]))) + 1,
                            value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
              }), br(), checkboxInput("times", "Predicted Survival at this Follow Up:"),
              conditionalPanel(condition = "input.times == true",
                               sliderInput(paste("numeric", (nn - 1), sep = ""), i.numeric[(nn - 1), 1],
                                           min = as.integer(min(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))),
                                           max = as.integer(max(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))) + 1,
                                           value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))))))
            }
            if (nn == 2){
              slide.bars <- list(br(), checkboxInput("times", "Predicted Survival at this Follow Up:"),
              conditionalPanel(condition = "input.times == true",
                               sliderInput(paste("numeric", (nn - 1), sep = ""), i.numeric[(nn - 1), 1],
                                           min = as.integer(min(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))),
                                           max = as.integer(max(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))) + 1,
                                           value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))))))
            }
          }

          if (covariate == "numeric") {
            if (nn > 2){
              slide.bars <- list(lapply(1:(nn - 2), function(j) {
                numericInput(paste("numeric", j, sep = ""), i.numeric[j, 1],
                             value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
              }), br(), checkboxInput("times", "Predicted Survival at this Follow Up:"),
              conditionalPanel(condition = "input.times == true",
                               numericInput(paste("numeric", (nn - 1), sep = ""), i.numeric[(nn - 1), 1],
                                            value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))))))
            }
            if (nn == 2){
              slide.bars <- list(br(), checkboxInput("times", "Predicted Survival at this Follow Up:"),
              conditionalPanel(condition = "input.times == true",
                               numericInput(paste("numeric", (nn - 1), sep = ""), i.numeric[(nn - 1), 1],
                                            value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))))))
            }
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
        if (nn > 1) {
          input.n <- vector("list", (nn - 1))
          for (i in 1:(nn - 1)) {
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
          out2 <- cbind(out[wher], NO=input$add)
          input.data <<- rbind(input.data[-1], out2)
        }
        if (a > 0) {
          wher <- match(names(out), names(input.data))
          out2 <- cbind(out[wher], NO=input$add)
          input.data <<- rbind(input.data, out2)
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
          OUT <- isolate({
            if (isTRUE(compare(old.d, new.d())) == FALSE) {
              new.d <- cbind(stat = 1, new.d())
              names(new.d)[1] <- paste(sts)
              if (n.strata > 0) {
                pred <- predict(model, newdata = new.d, se.fit = TRUE,
                                conf.int = clevel, type = "expected", reference = "strata")
              }
              if (n.strata == 0) {
                pred <- predict(model, newdata = new.d, se.fit = TRUE,
                                conf.int = clevel, type = "expected")
              }
              upb <- exp(-(pred$fit - (qnorm(1 - (1 - clevel)/2) * pred$se.fit)))
              if (upb > 1) {
                upb <- 1
              }
              lwb <- exp(-(pred$fit + (qnorm(1 - (1 - clevel)/2) * pred$se.fit)))
              if (ptype == "st") {
                d.p <- data.frame(Prediction = exp(-pred$fit), Lower.bound = lwb,
                                  Upper.bound = upb)
              }
              if (ptype == "1-st") {
                d.p <- data.frame(Prediction = 1 - exp(-pred$fit), Lower.bound = 1 - upb,
                                  Upper.bound = 1 - lwb)
              }

              old.d <<- new.d()
              data.p <- cbind(d.p, counter = 1, NO=input$add)
              p1 <<- rbind(p1, data.p)
              p1$count <- seq(1, dim(p1)[1])
              p1
            } else {
              p1$count <- seq(1, dim(p1)[1])
              OUT <- p1
            }
          })
        }
        OUT
      })

      s.fr <- NULL
      old.d2 <- NULL
      b <- 1
      St <- TRUE

      if (n.strata > 0) {
        sub.fit1 <- reactive({
          nam <- NULL
          aa <- 0
          fit1 <- survfit(model, newdata = new.d())
          l.s <- attr(model$terms, "specials")$strata
          for (i in l.s) {
            nam0 <- paste(new.d()[[which(i.factor[, 2] == i)]], sep = "")
            if (aa == 0) {
              nam <- paste(nam0)
            }
            if (aa > 0) {
              nam <- paste(nam, ", ", nam0, sep = "")
            }
            aa <- aa + 1
          }
          sub.fit1 <- subset(as.data.frame(summary(fit1)[2:8]), strata == nam)
          return(sub.fit1)
        })
      }

      dat.p <- reactive({
        if (isTRUE(compare(old.d2, new.d())) == FALSE) {
          s.frame <- isolate({
            fit1 <- survfit(model, newdata = new.d())
            if (n.strata == 0) {
              sff <- as.data.frame(summary(fit1)[2:8])
              sff <- cbind(sff, event=1-sff[[5]], part = b)
              if (sff$time[1] != 0){
                sff2 <- sff[1, ]
                sff2[1, ] <- NA
                sff2$time[1] <- 0
                sff2$n.risk[1] <- model$n
                sff2$surv[1] <- 1
                sff2$event[1] <- 0
                sff2$part[1] <- sff$part[1]
                s.f <- rbind(sff2, sff)
              } else {
                s.f <- sff
              }
            }
            if (n.strata > 0) {
              sff <- cbind(sub.fit1(), part = b)
              sff <- cbind(sff, event=1-sff[[6]])
              if (sff$time[1] != 0) {
                sff2 <- sff[1, ]
                sff2[1, ] <- NA
                sff2$time[1] <- 0
                sff2$n.risk[1] <- sff[1,2]
                sff2$surv[1] <- 1
                sff2$event[1] <- 0
                sff2$part[1] <- sff$part[1]
                s.f <- rbind(sff2, sff)
              } else {
                s.f <- sff
              }
              s.f$n.risk <- s.f$n.risk/s.f$n.risk[1]
            }
            if (dim(s.f)[1] < 3) {
              St <<- FALSE
              stop("Error in data structure: There is not enough data in the current strata level")
            }
            s.fr <<- rbind(s.fr, s.f)
            old.d2 <<- new.d()
            b <<- b + 1
            s.fr
          })
        } else {
          s.frame <- isolate({
            s.fr
          })
        }
      })

      output$plot <- renderPlot({
        if (St == TRUE) {
          if (input$add == 0)
            return(NULL)
          if (input$add > 0) {
            if (input$trans == TRUE) {
              if (ptype == "st") {
                pl <- isolate({
                  p2 <- ggplot(data = dat.p())
                  p2 <- p2 + geom_step(aes(x = time, y = surv, alpha = n.risk, color = as.factor(part), group = part))
                  p2 <- p2 + ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05)
                  p2 <- p2 + labs(title = "Estimated Survival Probability", x = "Follow Up Time", y = "S(t)") + theme_bw()
                  p2 <- p2 + theme(text = element_text(face = "bold", size = 14), legend.position = "none")
                })
              }
              if (ptype == "1-st") {
                pl <- isolate({
                  p2 <- ggplot(data = dat.p())
                  p2 <- p2 + geom_step(aes(x = time, y = event, alpha = n.risk, color = as.factor(part), group = part))
                  p2 <- p2 + ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05)
                  p2 <- p2 + labs(title = "Estimated Probability", x = "Follow Up Time", y = "F(t)")
                  p2 <- p2 + theme_bw() + theme(text = element_text(face = "bold", size = 14), legend.position = "none")
                })
              }
            }
            if (input$trans == FALSE) {
              if (ptype == "st") {
                pl <- isolate({
                  p2 <- ggplot(data = dat.p())
                  p2 <- p2 + geom_step(aes(x = time, y = surv, color = as.factor(part), group = part))
                  p2 <- p2 + ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05)
                  p2 <- p2 + labs(title = "Estimated Survival Probability", x = "Follow Up Time", y = "S(t)") + theme_bw()
                  p2 <- p2 + theme(text = element_text(face = "bold", size = 14), legend.position = "none")
                })
              }
              if (ptype == "1-st") {
                pl <- isolate({
                  p2 <- ggplot(data = dat.p())
                  p2 <- p2 + geom_step(aes(x = time, y = event, color = as.factor(part), group = part))
                  p2 <- p2 + ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05)
                  p2 <- p2 + labs(title = "Estimated Probability", x = "Follow Up Time", y = "F(t)")
                  p2 <- p2 + theme_bw() + theme(text = element_text(face = "bold", size = 14), legend.position = "none")
                })
              }
            }
          }
          data2()
          print(pl)
        }
        if (St == FALSE) {
          print("Restart the application")
        }
      })

      output$plot2 <- renderPlot({
        if (input$add == 0)
          return(NULL)
        isolate({
          if (is.null(new.d()))
            return(NULL)
          lim <- c(0, 1)
          yli <- c(0 - 0.5, 10 + 0.5)
          if (dim(input.data)[1] > 11)
            yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
          p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
          p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
          p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
          p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
                           size = 1.45, height = 0.4, colour = data2()$count)
          if (ptype == "st") {
            p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for Survival Probability", sep = ""),
                          x = "Survival Probability", y = NULL)
          }
          if (ptype == "1-st") {
            p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for F(t)", sep = ""),
                          x = "Probability", y = NULL)
          }
          p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 14))
          print(p)
        })
        data2()
      })

      output$data.pred <- renderPrint({
        if (input$add > 0) {
          isolate({
            if (nrow(data2() > 0)) {
              di <- ncol(input.data)
              data.p <- merge(input.data[-1, ], data2()[1:5], by="NO")
              data.p <- data.p[, !(colnames(data.p) %in% c("NO", "counter"))]
              stargazer(data.p, summary = FALSE, type = "text")
            }
          })
        }
      })

      output$summary <- renderPrint({
        coef.c <- exp(model$coef)
        ci.c <- exp(suppressMessages(confint(model, level = clevel)))
        stargazer(model, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F,
                  type = "text", omit.stat = c("LL", "ser", "f"), ci = TRUE, ci.level = clevel,
                  single.row = TRUE, title = paste("Cox model:", model$call[2], sep = " "))
      })
    }
  )
  )
}
