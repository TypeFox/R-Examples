data.Normalization<-function (x, type = "n0", normalization = "column") 
{
    bycolumn = T
    if (normalization == "row") 
        bycolumn = F
    if (is.vector(x) && !is.list(x)) {
        if (is.numeric(resul <- x)) 
            resul <- switch(type, n0 = x, n1 = (x - mean(x))/sd(x), 
                n2 = (x - median(x))/mad(x), n3 = (x - mean(x))/(max(x) - 
                  min(x)), n3a = (x - median(x))/(max(x) - min(x)), 
                n4 = (x - min(x))/(max(x) - min(x)), n5 = (x - 
                  mean(x))/(max(abs((x) - mean(x)))), n5a = (x - 
                  median(x))/(max(abs((x) - median(x)))), n6 = x/sd(x), 
                n6a = x/mad(x), n7 = x/(max(x) - min(x)), n8 = x/max(x), 
                n9 = x/mean(x), n9a = x/median(x), n10 = x/sum(x), 
                n11 = x/(sum(x^2)^0.5), n12 = (x - mean(x))/(sum((x - 
                  mean(x))^2)^0.5), 
                n12a = (x - median(x))/(sum((x - median(x))^2)^0.5),
                n13 = (x - ((max(x)+min(x))/2))/((max(x) - min(x))/2))
        else warning("Data not numeric, normalization not applicable")
        names(resul) <- names(x)
    }
    else if (is.data.frame(x)) {
        resul <- NULL
        if (bycolumn) {
            for (nn in names(x)) {
                if (is.numeric(x[, nn])) {
                  resul <- switch(type, n0 = cbind(resul, (x[, 
                    nn])), n1 = cbind(resul, (x[, nn] - mean(x[, 
                    nn]))/(sd(x[, nn]))), n2 = cbind(resul, (x[, 
                    nn] - median(x[, nn]))/(mad(x[, nn]))), n3 = cbind(resul, 
                    (x[, nn] - mean(x[, nn]))/(max(x[, nn]) - 
                      min(x[, nn]))), n3a = cbind(resul, (x[, 
                    nn] - median(x[, nn]))/(max(x[, nn]) - min(x[, 
                    nn]))), n4 = cbind(resul, (x[, nn] - min(x[, 
                    nn]))/(max(x[, nn]) - min(x[, nn]))), n5 = cbind(resul, 
                    (x[, nn] - mean(x[, nn]))/(max(abs(x[, nn] - 
                      mean(x[, nn]))))), n5a = cbind(resul, (x[, 
                    nn] - median(x[, nn]))/(max(abs(x[, nn] - 
                    median(x[, nn]))))), n6 = cbind(resul, (x[, 
                    nn])/sd(x[, nn])), n6a = cbind(resul, (x[, 
                    nn])/mad(x[, nn])), n7 = cbind(resul, (x[, 
                    nn])/(max(x[, nn]) - min(x[, nn]))), n8 = cbind(resul, 
                    (x[, nn])/(max(x[, nn]))), n9 = cbind(resul, 
                    (x[, nn])/(mean(x[, nn]))), n9a = cbind(resul, 
                    (x[, nn])/(median(x[, nn]))), n10 = cbind(resul, 
                    (x[, nn])/(sum(x[, nn]))), n11 = cbind(resul, 
                    (x[, nn])/(sum(x[, nn]^2)^0.5)), n12 = cbind(resul, 
                    (x[, nn] - mean(x[, nn]))/(sum((x[, nn] - 
                      mean(x[, nn]))^2)^0.5)), 
                     n12a = cbind(resul, (x[, nn] - median(x[, nn]))/(sum((x[, nn] - 
                      median(x[, nn]))^2)^0.5)),
                     n13 = cbind(resul, (x[, nn] - ((max(x[, nn])+min(x[, nn]))/2))/((max(x[, nn]) - min(x[, nn]))/2)))
                }
                else {
                  resul <- cbind(resul, x[, nn])
                  warning("Data not numeric, normalization not applicable")
                }
            }
        }
        else {
            for (nn in 1:nrow(x)) {
                if (sum(is.na(as.numeric((x[nn, ])))) == 0) {
                  resul <- switch(type, n0 = rbind(resul, (x[nn, 
                    ])), n1 = rbind(resul, (x[nn, ] - mean(as.numeric(x[nn, 
                    ])))/(sd(as.numeric(x[nn, ])))), n2 = rbind(resul, 
                    (x[nn, ] - median(as.numeric(x[nn, ])))/(mad(as.numeric(x[nn, 
                      ])))), n3 = rbind(resul, (x[nn, ] - mean(as.numeric(x[nn, 
                    ])))/(max(as.numeric(x[nn, ])) - min(as.numeric(x[nn, 
                    ])))), n3a = rbind(resul, (x[nn, ] - median(as.numeric(x[nn, 
                    ])))/(max(as.numeric(x[nn, ])) - min(as.numeric(x[nn, 
                    ])))), n4 = rbind(resul, (x[nn, ] - min(as.numeric(x[nn, 
                    ])))/(max(as.numeric(x[nn, ])) - min(as.numeric(x[nn, 
                    ])))), n5 = rbind(resul, (x[nn, ] - mean(as.numeric(x[nn, 
                    ])))/(max(abs(as.numeric(x[nn, ]) - mean(as.numeric(x[nn, 
                    ])))))), n5a = rbind(resul, (x[nn, ] - median(as.numeric(x[nn, 
                    ])))/(max(abs(as.numeric(x[nn, ]) - median(as.numeric(x[nn, 
                    ])))))), n6 = rbind(resul, (x[nn, ])/sd(as.numeric(x[nn, 
                    ]))), n6a = rbind(resul, (x[nn, ])/mad(as.numeric(x[nn, 
                    ]))), n7 = rbind(resul, (x[nn, ])/(max(as.numeric(x[nn, 
                    ])) - min(as.numeric(x[nn, ])))), n8 = rbind(resul, 
                    (x[nn, ])/(max(as.numeric(x[nn, ])))), n9 = rbind(resul, 
                    (x[nn, ])/(mean(as.numeric(x[nn, ])))), n9a = rbind(resul, 
                    (x[nn, ])/(median(as.numeric(x[nn, ])))), 
                    n10 = rbind(resul, (x[nn, ])/(sum(as.numeric(x[nn, 
                      ])))), n11 = rbind(resul, (x[nn, ])/(sum(as.numeric(x[nn, 
                      ])^2)^0.5)), n12 = rbind(resul, (x[nn, 
                      ] - mean(as.numeric(x[nn, ])))/(sum((as.numeric(x[nn, 
                      ]) - mean(as.numeric(x[nn, ])))^2)^0.5)), 
                    n12a = rbind(resul, (x[nn, ] - median(as.numeric(x[nn, 
                      ])))/(sum((as.numeric(x[nn, ]) - median(as.numeric(x[nn, 
                      ])))^2)^0.5)),
                    n13 = rbind(resul, (x[nn, ] - ((max(as.numeric(x[nn, ]))+min(as.numeric(x[nn, ])))/2))/
                    ((max(as.numeric(x[nn, ])) - min(as.numeric(x[nn, ])))/2)))
                }
                else {
                  resul <- rbind(resul, x[nn, ])
                  warning("Data not numeric, normalization not applicable")
                }
            }
        }
        resul <- data.frame(resul)
        names(resul) <- names(x)
        row.names(resul) <- row.names(x)
    }
    else if (is.matrix(x)) {
        if (is.numeric(resul <- x)) {
            resul <- NULL
            if (bycolumn) {
                for (i in 1:ncol(x)) resul <- switch(type, n0 = cbind(resul, 
                  (x[, i])), n1 = cbind(resul, (x[, i] - mean(x[, 
                  i]))/(sd(x[, i]))), n2 = cbind(resul, (x[, 
                  i] - median(x[, i]))/(mad(x[, i]))), n3 = cbind(resul, 
                  (x[, i] - mean(x[, i]))/(max(x[, i]) - min(x[, 
                    i]))), n3a = cbind(resul, (x[, i] - median(x[, 
                  i]))/(max(x[, i]) - min(x[, i]))), n4 = cbind(resul, 
                  (x[, i] - min(x[, i]))/(max(x[, i]) - min(x[, 
                    i]))), n5 = cbind(resul, (x[, i] - mean(x[, 
                  i]))/(max(abs(x[, i] - mean(x[, i]))))), n5a = cbind(resul, 
                  (x[, i] - median(x[, i]))/(max(abs(x[, i] - 
                    median(x[, i]))))), n6 = cbind(resul, (x[, 
                  i])/sd(x[, i])), n6a = cbind(resul, (x[, i])/mad(x[, 
                  i])), n7 = cbind(resul, (x[, i])/(max(x[, i]) - 
                  min(x[, i]))), n8 = cbind(resul, (x[, i])/(max(x[, 
                  i]))), n9 = cbind(resul, (x[, i])/(mean(x[, 
                  i]))), n9a = cbind(resul, (x[, i])/(median(x[, 
                  i]))), n10 = cbind(resul, (x[, i])/(sum(x[, 
                  i]))), n11 = cbind(resul, (x[, i])/(sum(x[, 
                  i]^2)^0.5)), n12 = cbind(resul, (x[, i] - mean(x[, 
                  i]))/(sum((x[, i] - mean(x[, i]))^2)^0.5)), 
                  n12a = cbind(resul, (x[, i] - median(x[, i]))/(sum((x[, 
                    i] - median(x[, i]))^2)^0.5)),
                  n13 = cbind(resul, (x[, i] - ((max(x[, i])+min(x[, i]))/2))/((max(x[, i]) - min(x[, i]))/2)))
            }
            else {
                for (i in 1:nrow(x)) resul <- switch(type, n0 = rbind(resul, 
                  (x[i, ])), n1 = rbind(resul, (x[i, ] - mean(x[i, 
                  ]))/(sd(x[i, ]))), n2 = rbind(resul, (x[i, 
                  ] - median(x[i, ]))/(mad(x[i, ]))), n3 = rbind(resul, 
                  (x[i, ] - mean(x[i, ]))/(max(x[i, ]) - min(x[i, 
                    ]))), n3a = rbind(resul, (x[i, ] - median(x[i, 
                  ]))/(max(x[i, ]) - min(x[i, ]))), n4 = rbind(resul, 
                  (x[i, ] - min(x[i, ]))/(max(x[i, ]) - min(x[i, 
                    ]))), n5 = rbind(resul, (x[i, ] - mean(x[i, 
                  ]))/(max(abs(x[i, ] - mean(x[i, ]))))), n5a = rbind(resul, 
                  (x[i, ] - median(x[i, ]))/(max(abs(x[i, ] - 
                    median(x[i, ]))))), n6 = rbind(resul, (x[i, 
                  ])/sd(x[i, ])), n6a = rbind(resul, (x[i, ])/mad(x[i, 
                  ])), n7 = rbind(resul, (x[i, ])/(max(x[i, ]) - 
                  min(x[i, ]))), n8 = rbind(resul, (x[i, ])/(max(x[i, 
                  ]))), n9 = rbind(resul, (x[i, ])/(mean(x[i, 
                  ]))), n9a = rbind(resul, (x[i, ])/(median(x[i, 
                  ]))), n10 = rbind(resul, (x[i, ])/(sum(x[i, 
                  ]))), n11 = rbind(resul, (x[i, ])/(sum(x[i, 
                  ]^2)^0.5)), n12 = rbind(resul, (x[i, ] - mean(x[i, 
                  ]))/(sum((x[i, ] - mean(x[i, ]))^2)^0.5)), 
                  n12a = rbind(resul, (x[i, ] - median(x[i, ]))/(sum((x[i, 
                    ] - median(x[i, ]))^2)^0.5)),
                  n13 = rbind(resul, (x[i, ] - ((max(x[i, ])+min(x[i, ]))/2))/((max(x[i, ]) - min(x[i, ]))/2)))
            }
        }
        else warning("Data not numeric, normalization not applicable")
        dimnames(resul) <- dimnames(x)
    }
    else if (is.list(x)) {
        resul <- list(length(x))
        xx <- as.numeric(x)
        for (i in 1:length(x)) if (is.numeric(resul[[i]] <- x[[i]])) 
            resul[[i]] <- switch(type, n0 = x[[i]], n1 = (x[[i]] - 
                mean(xx))/sd(xx), n2 = (x[[i]] - median(xx))/mad(xx), 
                n3 = (x[[i]] - mean(xx))/(max(xx) - min(xx)), 
                n3a = (x[[i]] - median(xx))/(max(xx) - min(xx)), 
                n4 = (x[[i]] - min(xx))/(max(xx) - min(xx)), 
                n5 = (x[[i]] - mean(xx))/(max(abs((xx) - mean(xx)))), 
                n5a = (x[[i]] - median(xx))/(max(abs((xx) - median(xx)))), 
                n6 = x[[i]]/sd(xx), n6a = x[[i]]/mad(xx), n7 = x[[i]]/(max(xx) - 
                  min(xx)), n8 = x[[i]]/(max(xx)), n9 = x[[i]]/(mean(xx)), 
                n9a = x[[i]]/(median(xx)), n10 = x[[i]]/(sum(xx)), 
                n11 = x[[i]]/(sum(xx^2)^0.5), n12 = (x[[i]] - 
                  mean(xx))/(sum((xx - mean(xx))^2)^0.5), n12a = (x[[i]] - 
                  median(xx))/(sum((xx - median(xx))^2)^0.5),
                n13 = (x[[i]] - ((max(xx)+min(xx))/2))/((max(xx) - min(xx))/2))
        else warning("Data not numeric, normalization not applicable")
    }
    else if (!is.numeric(resul <- x)) 
        warning("Data not numeric, normalization not applicable")
    else stop("unknown input type")
    resul
}
