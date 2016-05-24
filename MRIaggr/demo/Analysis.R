#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% test file 1.3 : analysis
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# setwd("V:/etude25/2891")

require(MRIaggr)

data("MRIaggr.Pat1_red", package="MRIaggr")
summary(MRIaggr.Pat1_red)
str(MRIaggr.Pat1_red,max.level = 2)
summary(MRIaggr.Pat1_red, history = TRUE, param = TRUE)

#### 0- pre-processing ####
res <- calcContralateral(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"), num = NULL, type = "1NN_penalised", param.ref = "T1_t0",
                         distband = 1, lambda = 1,
                         verbose = TRUE, overwrite = TRUE, update.object = TRUE)

res <- calcRegionalContrast(MRIaggr.Pat1_red, bandwidth = 1.875, param = c("DWI_t0","MTT_t0","TTP_t0"),
                       W.range = 6, W.spatial_res = c(1.875,1.875,6),
                       verbose = TRUE, overwrite = TRUE, update.object = TRUE)

res <- calcNormalization(MRIaggr.Pat1_red, param = c("DWI_t0","MTT_t0","TTP_t0","T2_FLAIR_t2",
                                                  "DWI_t0_regional","MTT_t0_regional","TTP_t0_regional"),
                         update.object = TRUE, overwrite = TRUE)

# multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2_regional")

#### 1- display the lesion ####

# 2D plot
multiplot(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2")

# 3D plot
plotLesion3D(MRIaggr.Pat1_red, mask = "MASK_T2_FLAIR_t2", spatial_res = c(1.875,1.875,6),
             numeric2logical = TRUE)

# Volume plot
plotTableLesion(MRIaggr.Pat1_red, num = 1:3, type = "evolution",
                 mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"))

plotTableLesion(MRIaggr.Pat1_red, num = 1:3, type = "matplot",
                 col = 1:10,
                 mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"))


#### 2- predictive value of contrast parameters ####

#### exploratory ####

# distribution according tissue class for DWI t1
plotDistClass(MRIaggr.Pat1_red, param = "DWI_t0", bw.adjust = 2,
              col = c("red","green","blue","black"),
              param.membership  = c("CSF","WM","GM","MASK_T2_FLAIR_t2"))

# boxplot for lesioned vs intact voxels for several parameters
boxplotMask(MRIaggr.Pat1_red, param = c("DWI_t0","TTP_t0","MTT_t0"), mask = "MASK_T2_FLAIR_t2",
            numeric2logical = TRUE)

# correlation map
heatmapMRIaggr(MRIaggr.Pat1_red, param = c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las = 2, type = "image.plot", cex = 0.5,
           breaks = seq(-1, 1, length.out = 51),
           col = cm.colors(50))

#### impact of normalisation ####
carto_test <- selectContrast(MRIaggr.Pat1_red, na.rm = TRUE,                          
                          param = c("DWI_t0","DWI_t0_contro","T2_FLAIR_t2","T2_FLAIR_t2_contro","MASK_DWI_t0"))

df.contrastMASK <- data.frame(matrix(NA, nrow = 3, ncol = 2))
colnames(df.contrastMASK) <- c("AUC","AUPRC")
rownames(df.contrastMASK) <- c("DWI_t0","DWI_t0_contro","DWI_t0andDWI_t0_contro")

model <- glm(MASK_DWI_t0 ~ DWI_t0, family = binomial(link = "logit"), data = carto_test)
df.contrastMASK["DWI_t0","AUPRC"] <- calcAUPRC(x = round(model$fitted, digit = 2), y = carto_test$MASK_DWI_t0)["AUPRC"]
if(require(pROC)){
df.contrastMASK["DWI_t0","AUC"] <- auc(roc(carto_test$MASK_DWI_t0~round(model$fitted, digit = 2)))
}

model <- glm(MASK_DWI_t0 ~ DWI_t0_contro, family=binomial(link="logit"),data = carto_test)
df.contrastMASK["DWI_t0_contro","AUPRC"] <- calcAUPRC(x = round(model$fitted, digit = 2), y = carto_test$MASK_DWI_t0)["AUPRC"]
if(require(pROC)){
  df.contrastMASK["DWI_t0_contro","AUC"] <- auc(roc(carto_test$MASK_DWI_t0 ~ round(model$fitted, digit = 2)))
}

model <- glm(MASK_DWI_t0 ~ DWI_t0 + DWI_t0_contro, family = binomial(link = "logit"),data = carto_test)
df.contrastMASK["DWI_t0andDWI_t0_contro","AUPRC"] <- calcAUPRC(x = round(model$fitted, digit = 2), y = carto_test$MASK_DWI_t0)["AUPRC"]
if(require(pROC)){
  df.contrastMASK["DWI_t0andDWI_t0_contro","AUC"] <- auc(roc(carto_test$MASK_DWI_t0 ~ round(model$fitted, digit = 2)))
}

df.contrastMASK

#### glm model  ####
param <- c( "DWI_t0","DWI_t0_regional","MTT_t0","TTP_t0","MTT_t0_regional","TTP_t0_regional")
data <- selectContrast(MRIaggr.Pat1_red, param = param, coords = TRUE,
                    hemisphere = "lesion",
                    norm_mu = "contralateral", norm_sigma = "contralateral")
data$MASK_T2_FLAIR_t2 <- as.logical(selectContrast(MRIaggr.Pat1_red,
                           param = "MASK_T2_FLAIR_t2", hemisphere = "lesion", format = "vector"))

glm_test <- glm(MASK_T2_FLAIR_t2 ~ DWI_t0 + TTP_t0 + MTT_t0+
                DWI_t0_regional + MTT_t0_regional + TTP_t0_regional,
                family = binomial(link = "logit"), data = data)

data$predictions <- predict(glm_test, type = "response")
data$residuals <- data$predictions - data$MASK_T2_FLAIR_t2

multiplot(data[,c("i","j","k")],
             contrast = data$predictions, legend = FALSE,
             main = "predictions - slice ", breaks = seq(0, 1, length.out = 100)
)

multiplot(data[,c("i","j","k")],
             contrast = data$residuals, legend = FALSE,
             palette = "cm.colors", main = "residuals - slice ",breaks = seq(-1, 1, length.out = 100)
)

