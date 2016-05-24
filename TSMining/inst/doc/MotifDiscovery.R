## ----_1, echo=FALSE------------------------------------------------------
library(TSMining)
data(BuildOperation)
summary(BuildOperation)

## ----_2------------------------------------------------------------------
res.wcc <- Func.motif(ts = BuildOperation$WCC, global.norm = T, local.norm = F, window.size = 24, overlap = 0, w = 6, a = 5, mask.size = 5, max.dist.ratio = 1.2, count.ratio.1 = 1.1, count.ratio.2 = 1.1)

res.ahu <- Func.motif(ts = BuildOperation$AHU, global.norm = T, local.norm = F, window.size = 24, overlap = 0, w = 6, a = 5, mask.size = 5, max.dist.ratio = 1.2, count.ratio.1 = 1.1, count.ratio.2 = 1.1)

## ----_3------------------------------------------------------------------
library(ggplot2)
#Visualization
data.wcc <- Func.visual.SingleMotif(single.ts = BuildOperation$WCC, window.size = 24, motif.indices = res.wcc$Indices)
data.ahu <- Func.visual.SingleMotif(single.ts = BuildOperation$AHU, window.size = 24, motif.indices = res.ahu$Indices)

#Determine the total number of motifs discovered in the time series of WCC
n <- length(unique(data.wcc$data.1$Y))
#Make the plot
ggplot(data = data.wcc$data.1) +  
    geom_line(aes(x = 1:dim(data.wcc$data.1)[1], y = X)) +
    geom_point(aes(x = 1:dim(data.wcc$data.1)[1], y = X, color=Y, shape=Y))+
    scale_shape_manual(values = seq(from = 1, to = n)) +
    guides(shape=guide_legend(nrow = 2)) +
    xlab("Time (15-min)") + ylab("WCC Power Consumption (kW)") +
    theme(panel.background=element_rect(fill = "white", colour = "black"),
          legend.position="top",
          legend.title=element_blank())

#Determine the total number of motifs discovered in the time series of AHU
n <- length(unique(data.ahu$data.1$Y))
#Make the plot
ggplot(data = data.ahu$data.1) +  
    geom_line(aes(x = 1:dim(data.ahu$data.1)[1], y = X)) +
    geom_point(aes(x = 1:dim(data.ahu$data.1)[1], y = X, color=Y, shape=Y))+
    scale_shape_manual(values = seq(from = 1, to = n)) +
    guides(shape=guide_legend(nrow = 2)) +
    xlab("Time (15-min)") + ylab("AHU Power Consumption (kW)") +
    theme(panel.background=element_rect(fill = "white", colour = "black"),
          legend.position="top",
          legend.title=element_blank())

## ----_4------------------------------------------------------------------
for(i in 1:length(data.wcc$data.2)) {
    data.temp <- data.wcc$data.2[[i]]
    print(ggplot(data = data.temp) +  
        geom_line(aes(x = Time, y = Value, color=Instance, linetype=Instance)) +
        xlab("Time (15-min)") + ylab("WCC Power Consumption (kW)") + ggtitle(paste0("WCC Motif ",i)) +
        scale_y_continuous(limits=c(0,max(data.temp$Value))) +
        theme(panel.background=element_rect(fill = "white", colour = "black"),
              legend.position="none",
              legend.title=element_blank()))    
}

## ----_5------------------------------------------------------------------
for(i in 1:length(data.ahu$data.2)) {
    data.temp <- data.ahu$data.2[[i]]
    print(ggplot(data = data.temp) +  
              geom_line(aes(x = Time, y = Value, color=Instance, linetype=Instance)) +
              xlab("Time (15-min)") + ylab("AHU Power Consumption (kW)") + ggtitle(paste0("AHU Motif ",i)) +
              scale_y_continuous(limits=c(0,max(data.temp$Value))) +
              theme(panel.background=element_rect(fill = "white", colour = "black"),
                    legend.position="none",
                    legend.title=element_blank()))    
}

## ----_6------------------------------------------------------------------
res.multi <- Func.motif.multivariate(motif.list = list(res.wcc$Indices, res.ahu$Indices), window.sizes = c(24,24), alpha = .7)

## ----_7------------------------------------------------------------------
#Focus on the third multivariate motif
data.multi <- Func.visual.MultiMotif(data = BuildOperation[,c("WCC","AHU")], multi.motifs = res.multi, index = 3)

ggplot(data = data.multi, aes(x = T, y = X)) + geom_line() + geom_point(aes(col=Lab, shape=Lab)) + facet_grid(Facet~.) +
    theme(panel.background=element_rect(fill = "white", colour = "black"), 
          legend.title=element_blank(),
          legend.position="top")

