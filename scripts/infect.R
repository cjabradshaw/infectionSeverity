## ICU infections

## libraries
library(boot)
library(dismo)
library(dplyr)
library(gbm)
library(ggplot2)
library(gridExtra)
library(usdm)

# custom functions
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

modifyVecFunc <- function(obj_name, index, new_value) {
  tryCatch({
    if (!object_exists(obj_name)) {
      stop("object does not exist: ", obj_name)
    }
    
    obj <- get(obj_name)
    
    # check object type and handle accordingly
    if (is.vector(obj) && !is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[index] <- new_value
    } else if (is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[[index]] <- new_value
    } else if (is.data.frame(obj)) {
      if (!all(index <= dim(obj))) stop("index out of bounds")
      obj[index[1], index[2]] <- new_value
    } else {
      stop("unsupported object type")
    }
    
    # Save modified object back to its original name
    assign(obj_name, obj, envir = .GlobalEnv)
    
  }, error = function(e) {
    message("error: ", e$message)
    return(FALSE)
  })
}

object_exists <- function(obj_name) {
  exists(obj_name, envir = .GlobalEnv)
}

## import data
infect.dat <- read.csv("infect.csv", header = T)
head(infect.dat)

## check columns
table(infect.dat$infect_grp)
hist(infect.dat$age)
table(infect.dat$gender)
hist(infect.dat$height)
hist(infect.dat$weight)
hist(infect.dat$AP3)
hist(infect.dat$GCS)
table(infect.dat$diabetes)
table(infect.dat$frailty)
hist(infect.dat$wcc_max)
hist(infect.dat$wcc_min)
hist(infect.dat$albumin)
hist(infect.dat$globulin)

## create response variables
infect.dat$infect_bin1 <- (ifelse(infect.dat$infect_grp == "CTRL", 0,
                                  ifelse(infect.dat$infect_grp == "OTHER", NA, 1)))
table(infect.dat$infect_bin1)
infect.dat$infect_bin2 <- (ifelse(infect.dat$infect_grp == "CTRL" | infect.dat$infect_grp == "OTHER", 0, 1))
table(infect.dat$infect_bin2)

## recode infection group to 3-letter code
infect.dat$infect_grp_code <- ifelse(infect.dat$infect_grp == "CTRL", "CTL",
                                      ifelse(infect.dat$infect_grp == "OTHER", "OTH",
                                             ifelse(infect.dat$infect_grp == "ASP PNEUMO", "ASP",
                                                    ifelse(infect.dat$infect_grp == "SEPSIS", "SEP",
                                                           ifelse(infect.dat$infect_grp == "VIRAL PNEUMO", "VIR", "BAC")))))
table(infect.dat$infect_grp_code)

## factorise gender
infect.dat$gender_fac <- factor(infect.dat$gender, levels = c("Male", "Female"))
infect.dat$gender_fac <- ifelse(infect.dat$gender_fac == "Male", "male",
                                 ifelse(infect.dat$gender_fac == "Female", "female", NA))
infect.dat$gender_fac <- factor(infect.dat$gender_fac, levels = c("male", "female"))
table(infect.dat$gender_fac)

## create body mass index
infect.dat$bmi <- infect.dat$weight / (infect.dat$height/100)^2
hist(log10(infect.dat$bmi))

## create frailty ordinal factor
infect.dat$frailty_ord <- factor(infect.dat$frailty,
                                 levels = c("VeryFit", "Well","ManagingWell","Vulnerable","MildlyFrail",
                                            "ModeratelyFrail","SeverelyFrail","ExtremelyFrail"), ordered = TRUE)
table(infect.dat$frailty_ord)

## create numeric ordinal frailty variable
infect.dat$frailty_num <- factor(as.numeric(infect.dat$frailty_ord), levels=1:8, ordered=TRUE)
table(infect.dat$frailty_num)
str(infect.dat)

## create frailty binary
infect.dat$frailty_bin <- as.factor(ifelse(infect.dat$frailty %in% c("VeryFit", "Well","ManagingWell"), 0, 1))
table(infect.dat$frailty_bin)

## create violin plot of infection group by age
ggplot(infect.dat, aes(x=infect_grp_code, y=age)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width=0.1, fill="white") +
  xlab("infection group") +
  ylab("age") +
  theme_minimal()

## violin plot of infection group by AP3
ggplot(infect.dat, aes(x=infect_grp_code, y=AP3)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width=0.1, fill="white") +
  xlab("infection group") +
  ylab("AP3 score") +
  theme_minimal()

## violin plot of infection group by GCS
ggplot(infect.dat, aes(x=infect_grp_code, y=GCS)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width=0.1, fill="white") +
  xlab("infection group") +
  ylab("GCS score") +
  theme_minimal()

## violin plot of infection group by wcc_med
which(is.na(infect.dat$wcc_max)==T)
which(is.na(infect.dat$wcc_min)==T)
infect.dat$wcc_med <- (infect.dat$wcc_max + infect.dat$wcc_min) / 2
ggplot(infect.dat, aes(x=infect_grp_code, y=wcc_med)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width=0.1, fill="white") +
  xlab("infection group") +
  ylab("median white cell count") +
  theme_minimal()
  
## violin plot of AP3 score vs. frailty_num
ggplot(infect.dat, aes(x=frailty_num, y=AP3)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width=0.1, fill="white") +
  xlab("frailty score (1=fit ... 8=extremely frail)") +
  ylab("AP3 score") +
  theme_minimal()

head(infect.dat)
str(infect.dat)

## logarithm of bmi
infect.dat$lbmi <- log10(infect.dat$bmi)
hist(infect.dat$lbmi)

## detect bmi outliers (modified z-score)
bmi_NA_sub <- which(is.na(infect.dat$lbmi)==T)
infect.dat.noNA <- infect.dat[-bmi_NA_sub, ]
med_lbmi <- median(infect.dat.noNA$lbmi, na.rm=T)
mad_lbmi <- median(abs(infect.dat.noNA$lbmi - med_lbmi), na.rm=T)
mzs_lbmi <- 0.6745 * (infect.dat.noNA$lbmi - med_lbmi) / mad_lbmi
lbmi_outliers <- which(abs(mzs_lbmi) > 3.5)

## remove bmi outliers
infect.dat2 <- infect.dat.noNA[-lbmi_outliers, ]
dim(infect.dat)[1] - dim(infect.dat.noNA)[1]
dim(infect.dat.noNA)[1] - dim(infect.dat2)[1]

## histogram of outlier-culled log bmi
hist(infect.dat2$lbmi)

# create log-transformed wcc_med variable
infect.dat2$lwcc_med <- log10(infect.dat2$wcc_med)

## examine globulin and albumin distributions
hist(infect.dat2$globulin)
hist(infect.dat2$albumin)

## scale  variables
infect.dat2$age.sc <- scale(infect.dat2$age, center = TRUE, scale = TRUE)
hist(infect.dat2$age.sc)
infect.dat2$lbmi.sc <- scale(infect.dat2$lbmi, center = TRUE, scale = TRUE)
hist(infect.dat2$lbmi.sc)
infect.dat2$AP3.sc <- scale(infect.dat2$AP3, center = TRUE, scale = TRUE)
hist(infect.dat2$AP3.sc)
infect.dat2$lwcc_med.sc <- scale(infect.dat2$lwcc_med, center = TRUE, scale = TRUE)
hist(infect.dat2$lwcc_med.sc)
infect.dat2$albumin.sc <- scale(infect.dat2$albumin, center = TRUE, scale = TRUE)
infect.dat2$globulin.sc <- scale(infect.dat2$globulin, center = TRUE, scale = TRUE)
str(infect.dat2)

## correlation matrix
correlation_vars <- c("age.sc","lbmi.sc","AP3.sc","lwcc_med.sc","albumin.sc","globulin.sc")
correlation_cols <- which(colnames(infect.dat2) %in% correlation_vars)
cor.mat <- cor(na.omit(infect.dat2[, correlation_vars]), use="pairwise.complete.obs", method="spearman")
cor.mat[upper.tri(cor.mat, diag=T)] <- NA
cor.mat

## variance inflation factor
vif(na.omit(infect.dat2[, correlation_vars]))

## bivariate plots
colnames(infect.dat2)
ggplot(infect.dat2, aes(x=AP3.sc, y=GCS)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  xlab("acute physiol, age, chronic health score") +
  ylab("Glasgow coma scale") +
  theme_minimal()

ggplot(infect.dat2, aes(x=AP3.sc, y=lwcc_med.sc)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  xlab("acute physiol, age, chronic health score") +
  ylab("white blood cell count") +
  theme_minimal()

ggplot(infect.dat2, aes(x=AP3.sc, y=albumin.sc)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  xlab("acute physiol, age, chronic health score") +
  ylab("blood albumin") +
  theme_minimal()

ggplot(infect.dat2, aes(x=AP3.sc, y=globulin.sc)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  xlab("acute physiol, age, chronic health score") +
  ylab("blood globulin") +
  theme_minimal()

ggplot(infect.dat2, aes(x=albumin.sc, y=globulin.sc)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color="red") +
  xlab("blood albumin") +
  ylab("blood globulin") +
  theme_minimal()

## response & predictors
## probability of infection ('other' included in non-infected group)
response_var <- "infect_bin2"
response_col <- which(colnames(infect.dat2) == response_var)
predictor_vars <- c("gender_fac","age.sc","lbmi.sc","frailty_num")
predictor_cols <- which(colnames(infect.dat2) %in% predictor_vars)
weights_var <- c("lwcc_med.sc") # lwcc_med.sc, AP3.sc, globulin.sc, albumin.sc
weights_col <- which(colnames(infect.dat2) %in% weights_var)

## create no-NA dataset with response, predictors, and weights for boosted regression tree
infect.dat2.brt <- infect.dat2[, c(response_col, predictor_cols, weights_col)]
infect.dat2.brt <- na.omit(infect.dat2.brt)
table(infect.dat2.brt$infect_bin2)
str(infect.dat2.brt)
head(infect.dat2.brt)
dim(infect.dat2.brt)

## new column indicators
predictor_cols_brt <- which(colnames(infect.dat2.brt) %in% predictor_vars)
response_col_brt <- which(colnames(infect.dat2.brt) == response_var)
weights_col_brt <- which(colnames(infect.dat2.brt) %in% weights_var)

## boosted regression tree
brt.bin2 <- gbm.step(infect.dat2.brt, gbm.x = attr(infect.dat2.brt, "names")[predictor_cols_brt],
                          gbm.y = attr(infect.dat2.brt, "names")[response_col_brt], 
                          weights=attr(infect.dat2.brt, "names")[weights_col_brt],
                          family="bernoulli", max.trees=100000,
                          tolerance = 0.00004, learning.rate = 0.00008, bag.fraction=0.75,
                          tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(brt.bin2)
barplot(summary(brt.bin2)$rel.inf, names.arg = summary(brt.bin2)$var, xlab="relative influence", ylab="", col="blue")
brt.bin2.summ <- summary(brt.bin2)

brt.bin2.CV.cor <- 100 * brt.bin2$cv.statistics$correlation.mean
brt.bin2.CV.cor.se <- 100 * brt.bin2$cv.statistics$correlation.se
print(c(brt.bin2.CV.cor, brt.bin2.CV.cor.se))

gbm.plot.fits(brt.bin2, v=0)

gbm.plot(brt.bin2, variable.no=0, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", plot.layout=c(2,2))
brt.bin2$gbm.call$gbm.x


gbm.plot(brt.bin2, variable.no=4, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="scaled log10 BMI", plot.layout=c(1,1))
gbm.plot(brt.bin2, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="scaled age (years)", plot.layout=c(1,1))
gbm.plot(brt.bin2, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="frailty score (1=fit ... 8=extremely frail)", plot.layout=c(1,1))
gbm.plot(brt.bin2, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="patient gender", plot.layout=c(1,1))

# bmi
bmi.pred.out <- plot.gbm(brt.bin2, i.var=4, return.grid=T)
bmi.pred.out$lbmi <- bmi.pred.out$lbmi.sc * attr(infect.dat2$lbmi.sc, "scaled:scale") + attr(infect.dat2$lbmi.sc, "scaled:center")
bmi.pred.out$bmi <- 10^(bmi.pred.out$lbmi)
bmi.pred.out$prob_inf <- logit2prob(bmi.pred.out$y)
head(bmi.pred.out)

# age
age.pred.out <- plot.gbm(brt.bin2, i.var=3, return.grid=T)
age.pred.out$age <- age.pred.out$age.sc * attr(infect.dat2$age.sc, "scaled:scale") + attr(infect.dat2$age.sc, "scaled:center")
age.pred.out$prob_inf <- logit2prob(age.pred.out$y)
head(age.pred.out)

## frailty
frailty.pred.out <- plot.gbm(brt.bin2, i.var=2, return.grid=T)
frailty.pred.out$prob_inf <- logit2prob(frailty.pred.out$y)
frailty.pred.out$frailty_ord <- levels(infect.dat2$frailty_ord)
frailty.pred.out

# gender
gender.pred.out <- plot.gbm(brt.bin2, i.var=1, return.grid = T)
gender.pred.out$prob_inf <- logit2prob(gender.pred.out$y)
gender.pred.out$gender <- levels(infect.dat2$gender_fac)
gender.pred.out

max.prob_inf <- 1.05 * max(c(bmi.pred.out$prob_inf, age.pred.out$prob_inf,
                      frailty.pred.out$prob_inf, gender.pred.out$prob_inf), na.rm=T)
min.prob_inf <- 0.95 * min(c(bmi.pred.out$prob_inf, age.pred.out$prob_inf,
                      frailty.pred.out$prob_inf, gender.pred.out$prob_inf), na.rm=T)

## 2 x 2 plot
## add predicted probabilities to original data for plotting
pred_pr_inf <- predict(brt.bin2, newdata=infect.dat2, type="response") 
infect.dat2$pred_pr_inf <- pred_pr_inf
hist(infect.dat2$pred_pr_inf)

plt.bmi <- ggplot(bmi.pred.out, aes(x=bmi, y=prob_inf)) +
  geom_point(data=infect.dat2, aes(x=bmi, y=pred_pr_inf), color="darkgrey", alpha=0.5) +
  geom_line(color="blue", linewidth=1.9) +
  xlab(expression("patient body mass index (kg/m"^2*")")) +
  ylab("predicted probability of infection") +
  ylim(min.prob_inf, max.prob_inf) +
  theme_minimal()

plt.age <- ggplot(age.pred.out, aes(x=age, y=prob_inf)) +
  geom_point(data=infect.dat2, aes(x=age, y=pred_pr_inf), color="darkgrey", alpha=0.5) +
  geom_line(color="blue", linewidth=1.9) +
  xlab("patient age (years)") +
  ylab("") +
  ylim(min.prob_inf, max.prob_inf) +
  theme_minimal()

frailty_dat <- na.omit(data.frame(frailty_num=infect.dat2$frailty_num, pred_pr_inf=infect.dat2$pred_pr_inf))
plt.frailty <- ggplot(frailty.pred.out, aes(x=frailty_num, y=prob_inf)) +
  geom_bar(stat="identity", fill="lightblue") +
  geom_violin(data=frailty_dat, aes(x=frailty_num, y=pred_pr_inf), fill="darkgrey", alpha=0.6) +
  scale_y_continuous(limits=c(0, max.prob_inf)) +
  xlab("patient admission frailty (1=fit ... 8=extremely frail)") +
  ylab("predicted probability of infection") +
  theme_minimal()

gender_dat <- na.omit(data.frame(gender_fac=infect.dat2$gender_fac, pred_pr_inf=infect.dat2$pred_pr_inf))
plt.gender <- ggplot(gender.pred.out, aes(x=gender_fac, y=prob_inf)) +
  geom_bar(stat="identity", fill="lightblue") +
  geom_violin(data=gender_dat, aes(x=gender_fac, y=pred_pr_inf), fill="darkgrey", alpha=0.6) +
  scale_y_continuous(limits=c(0, max.prob_inf)) +
  xlab("patient gender") +
  ylab("") +
  theme_minimal()

grid.arrange(plt.bmi, plt.age, plt.frailty, plt.gender, nrow=2, ncol=2)


## redo with 'other' included in infected group (infect_bin1)
# remove infect_bin1 values == NA
## create no-NA dataset with response, predictors, and weights for boosted regression tree
response_var1 <- "infect_bin1"
response_col1 <- which(colnames(infect.dat2) == response_var1)

infect.dat1.brt <- infect.dat2[, c(response_col1, predictor_cols, weights_col)]
infect.dat1.brt <- infect.dat1.brt[-which(is.na(infect.dat1.brt$infect_bin1)==T), ]
table(infect.dat1.brt$infect_bin1)
head(infect.dat1.brt)
str(infect.dat1.brt)
dim(infect.dat1.brt)

## new column indicators
response_col_brt1 <- which(colnames(infect.dat1.brt) == response_var1)
predictor_cols_brt1 <- which(colnames(infect.dat1.brt) %in% predictor_vars)
weights_col_brt1 <- which(colnames(infect.dat1.brt) %in% weights_var) 

## boosted regression tree (remove weights because too few data otherwise)
brt.bin1 <- gbm.step(infect.dat1.brt, gbm.x = attr(infect.dat1.brt, "names")[predictor_cols_brt1],
                     gbm.y = attr(infect.dat1.brt, "names")[response_col_brt1], 
                     #weights=attr(infect.dat1.brt, "names")[weights_col_brt1],
                     family="bernoulli", max.trees=100000,
                     tolerance = 0.00006, learning.rate = 0.00003, bag.fraction=0.75,
                     tree.complexity = 2, silent=F, tolerance.method = "auto")
summary(brt.bin1)
barplot(summary(brt.bin1)$rel.inf, names.arg = summary(brt.bin1)$var, xlab="relative influence", ylab="", col="blue")
brt.bin1.summ <- summary(brt.bin1)

brt.bin1.CV.cor <- 100 * brt.bin1$cv.statistics$correlation.mean
brt.bin1.CV.cor.se <- 100 * brt.bin1$cv.statistics$correlation.se
print(c(brt.bin1.CV.cor, brt.bin1.CV.cor.se))

gbm.plot(brt.bin1, variable.no=0, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", plot.layout=c(2,2))
brt.bin1$gbm.call$gbm.x

gbm.plot(brt.bin1, variable.no=4, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="scaled log10 BMI", plot.layout=c(1,1))
gbm.plot(brt.bin1, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="scaled age (years)", plot.layout=c(1,1))
gbm.plot(brt.bin1, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="frailty score (1=fit ... 8=extremely frail)", plot.layout=c(1,1))
gbm.plot(brt.bin1, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="Pr(infection)", x.label="patient gender", plot.layout=c(1,1))

# bmi
bmi.pred.out1 <- plot.gbm(brt.bin1, i.var=4, return.grid=T)
bmi.pred.out1$lbmi <- bmi.pred.out1$lbmi.sc * attr(infect.dat2$lbmi.sc, "scaled:scale") + attr(infect.dat2$lbmi.sc, "scaled:center")
bmi.pred.out1$bmi <- 10^(bmi.pred.out1$lbmi)
bmi.pred.out1$prob_inf <- logit2prob(bmi.pred.out1$y)
head(bmi.pred.out1)

# age
age.pred.out1 <- plot.gbm(brt.bin1, i.var=3, return.grid=T)
age.pred.out1$age <- age.pred.out1$age.sc * attr(infect.dat2$age.sc, "scaled:scale") + attr(infect.dat2$age.sc, "scaled:center")
age.pred.out1$prob_inf <- logit2prob(age.pred.out1$y)
head(age.pred.out1)

## frailty
frailty.pred.out1 <- plot.gbm(brt.bin1, i.var=2, return.grid=T)
frailty.pred.out1$prob_inf <- logit2prob(frailty.pred.out1$y)
frailty.pred.out1$frailty_ord <- levels(infect.dat2$frailty_ord)
frailty.pred.out1

# gender
gender.pred.out1 <- plot.gbm(brt.bin1, i.var=1, return.grid = T)
gender.pred.out1$prob_inf <- logit2prob(gender.pred.out1$y)
gender.pred.out1$gender <- levels(infect.dat2$gender_fac)
gender.pred.out1

max.prob_inf1 <- max(c(bmi.pred.out1$prob_inf, age.pred.out1$prob_inf,
                      frailty.pred.out1$prob_inf, gender.pred.out1$prob_inf), na.rm=T)
min.prob_inf1 <- min(c(bmi.pred.out1$prob_inf, age.pred.out1$prob_inf,
                      frailty.pred.out1$prob_inf, gender.pred.out1$prob_inf), na.rm=T)

## 2 x 2 plot
## add predicted probabilities to original data for plotting
pred_pr_inf1 <- predict(brt.bin1, newdata=infect.dat2, type="response") 
infect.dat2$pred_pr_inf1 <- pred_pr_inf1
hist(infect.dat2$pred_pr_inf1)

plt.bmi1 <- ggplot(bmi.pred.out1, aes(x=bmi, y=prob_inf)) +
  geom_point(data=infect.dat2, aes(x=bmi, y=pred_pr_inf1), color="darkgrey", alpha=0.5) +
  geom_line(color="blue", linewidth=1.9) +
  xlab(expression("patient body mass index (kg/m"^2*")")) +
  ylab("predicted probability of infection") +
  ylim(min.prob_inf1, 0.95) +
  theme_minimal()

plt.age1 <- ggplot(age.pred.out1, aes(x=age, y=prob_inf)) +
  geom_point(data=infect.dat2, aes(x=age, y=pred_pr_inf1), color="darkgrey", alpha=0.5) +
  geom_line(color="blue", linewidth=1.9) +
  xlab("patient age (years)") +
  ylab("") +
  ylim(min.prob_inf1, 0.95) +
  theme_minimal()

frailty_dat1 <- na.omit(data.frame(frailty_num=infect.dat2$frailty_num, pred_pr_inf1=infect.dat2$pred_pr_inf1))
plt.frailty1 <- ggplot(frailty.pred.out1, aes(x=frailty_num, y=prob_inf)) +
  geom_bar(stat="identity", fill="lightblue") +
  geom_violin(data=frailty_dat1, aes(x=frailty_num, y=pred_pr_inf1), fill="darkgrey", alpha=0.6) +
  scale_y_continuous(limits=c(0, 0.95)) +
  xlab("patient admission frailty (1=fit ... 8=extremely frail)") +
  ylab("predicted probability of infection") +
  theme_minimal()

gender_dat1 <- na.omit(data.frame(gender_fac=infect.dat2$gender_fac, pred_pr_inf1=infect.dat2$pred_pr_inf1))
plt.gender1 <- ggplot(gender.pred.out1, aes(x=gender_fac, y=prob_inf)) +
  geom_bar(stat="identity", fill="lightblue") +
  geom_violin(data=gender_dat1, aes(x=gender_fac, y=pred_pr_inf1), fill="darkgrey", alpha=0.6) +
  scale_y_continuous(limits=c(0, 0.95)) +
  xlab("patient gender") +
  ylab("") +
  theme_minimal()

grid.arrange(plt.bmi1, plt.age1, plt.frailty1, plt.gender1, nrow=2, ncol=2)



##############################################################################################
## resampled BRT for uncertainty in the infection severity (white blood cell count) weighting
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()
eq.sp.pts <- 100

head(infect.dat2)

## reset response & predictor cols
response_var <- "infect_bin2"
response_col <- which(colnames(infect.dat2) == response_var)
predictor_vars <- c("gender_fac","age.sc","lbmi.sc","frailty_num")
predictor_cols <- which(colnames(infect.dat2) %in% predictor_vars)

## create base dataset for resampling
wcc_cols <- which(colnames(infect.dat2) %in% c("wcc_min", "wcc_max"))
dat.rsmp <- na.omit(infect.dat2[, c(response_col, predictor_cols, wcc_cols)])
head(dat.rsmp)
dim(dat.rsmp)
table(dat.rsmp$infect_bin2)
dat.rsmp.orig <- dat.rsmp

## reset response & predictor cols for resampled dataset
response_col <- which(colnames(dat.rsmp) == response_var)
predictor_vars <- c("gender_fac","age.sc","lbmi.sc","frailty_num")
predictor_cols <- which(colnames(dat.rsmp) %in% predictor_vars)
npredictcols <- length(predictor_cols)

cont_predictors <- c("age.sc","lbmi.sc")
cont_predictor_cols <- which(brt.it$gbm.call$gbm.x %in% cont_predictors)
cont.seq <- seq(cont_predictor_cols[1], cont_predictor_cols[2], 1)
cont.it <- seq(1, length(cont_predictor_cols), 1)
iter.subtract <- length(predictor_cols) - length(cont_predictor_cols)

fac_predictor1 <- c("gender_fac")
fac_predictor1_col <- which(brt.it$gbm.call$gbm.x %in% fac_predictor1)

fac_predictor2 <- c("frailty_num")
fac_predictor2_col <- which(brt.it$gbm.call$gbm.x %in% fac_predictor2)

# create storage arrays
cont.val.arr <- cont.pred.arr <- array(data=NA, dim=c(eq.sp.pts, length(cont_predictors), biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), cont_predictors, paste("b",1:biter,sep="")))
fac1.pred.mat <- matrix(data=NA, nrow=biter, ncol=length(levels(infect.dat2$gender_fac)),
                             dimnames=list(paste("b",1:biter,sep=""), levels(infect.dat2$gender_fac)))
fac2.pred.mat <- matrix(data=NA, nrow=biter, ncol=length(levels(infect.dat2$frailty_num)),
                             dimnames=list(paste("b",1:biter,sep=""), levels(infect.dat2$frailty_num)))

# create storage vectors
ri.vec.names <- paste(predictor_vars,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:npredictcols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # reset dat.rsmp
  dat.rsmp <- dat.rsmp.orig
  
  ## randomly resample the weights (white blood cell count) for each patient
  dat.rsmp$wcc_ran <- runif(nrow(dat.rsmp), min=dat.rsmp$wcc_min, max=dat.rsmp$wcc_max)
  head(dat.rsmp)
  
  ## transform resampled wcc variable to log10 scale and scale for weighting
  dat.rsmp$lwcc_ran <- log10(dat.rsmp$wcc_ran)
  dat.rsmp$lwcc_ran.sc <- scale(dat.rsmp$lwcc_ran, center = TRUE, scale = TRUE)
  
  ## choose weighting column
  weights_var <- c("lwcc_ran.sc")
  weights_col <- which(colnames(dat.rsmp) == weights_var)
  
  ## boosted regression tree
  brt.it <- gbm.step(dat.rsmp, gbm.x = attr(dat.rsmp, "names")[predictor_cols],
                       gbm.y = attr(dat.rsmp, "names")[response_col], 
                       weights=attr(dat.rsmp, "names")[weights_col],
                       family="bernoulli", max.trees=100000,
                       tolerance = 0.0004, learning.rate = 0.0001, bag.fraction=0.75,
                       tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.it)==F) {
    brt.it.old <- brt.it
  }
  if (is.null(brt.it) == T) {
    brt.it <- gbm.step(dat.rsmp, gbm.x = attr(dat.rsmp, "names")[predictor_cols],
                       gbm.y = attr(dat.rsmp, "names")[response_col], 
                       weights=attr(dat.rsmp, "names")[weights_col],
                       family="bernoulli", max.trees=100000,
                       tolerance = 0.00004, learning.rate = 0.00008, bag.fraction=0.75,
                       tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.it) == T) {
    brt.it <- brt.it.old
  }
  
  # summary
  summ.fit <- summary(brt.it)
  
  if (is.null(brt.it) == F) {
    brt.it.old <- brt.it
  }
  
  # variable relative importance
  for (ri in 1:npredictcols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == predictor_vars[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.it$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.it$cv.statistics$correlation.se
  
  # response curves for continuous predictor variables (age, bmi)
  RESP.cont.val <- RESP.cont.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=length(cont_predictor_cols))
  for (p in cont.seq) {
    RESP.cont.val[,p-iter.subtract] <- plot.gbm(brt.it, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.cont.pred[,p-iter.subtract] <- plot.gbm(brt.it, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.cont.val.dat <- as.data.frame(RESP.cont.val)
  colnames(RESP.cont.val.dat) <- cont_predictors
  RESP.cont.pred.dat <- as.data.frame(RESP.cont.pred)
  colnames(RESP.cont.pred.dat) <- cont_predictors
  
  # response curves for factor predictor variables (gender_fac, frailty_num)
  # gender
  RESP.fac1.val <- plot.gbm(brt.it, i.var=fac_predictor1_col, return.grid=T)[,1]
  RESP.fac1.pred <- plot.gbm(brt.it, i.var=fac_predictor1_col, return.grid=T)[,2]
  RESP.fac1.val.dat <- as.data.frame(RESP.fac1.val)
  colnames(RESP.fac1.val.dat) <- fac_predictor1
  RESP.fac1.pred.dat <- as.data.frame(RESP.fac1.pred)
  colnames(RESP.fac1.pred.dat) <- fac_predictor1
  rownames(RESP.fac1.pred.dat) <- levels(infect.dat2$gender_fac)
  
  # frailty
  RESP.fac2.val <- plot.gbm(brt.it, i.var=fac_predictor2_col, return.grid=T)[,1]
  RESP.fac2.pred <- plot.gbm(brt.it, i.var=fac_predictor2_col, return.grid=T)[,2]
  RESP.fac2.val.dat <- as.data.frame(RESP.fac2.val)
  colnames(RESP.fac2.val.dat) <- fac_predictor2
  RESP.fac2.pred.dat <- as.data.frame(RESP.fac2.pred)
  colnames(RESP.fac2.pred.dat) <- fac_predictor2
  rownames(RESP.fac2.pred.dat) <- levels(infect.dat2$frailty_num)
                                         
  # add to storage arrays
  cont.val.arr[, , b] <- as.matrix(RESP.cont.val.dat)
  cont.pred.arr[, , b] <- as.matrix(RESP.cont.pred.dat)
  fac1.pred.mat[b, ] <- as.matrix(RESP.fac1.pred.dat)
  fac2.pred.mat[b, ] <- as.matrix(RESP.fac2.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2

# continuous variables
kappa.cont.n <- length(cont_predictors)
pred.cont.update <- cont.pred.arr[,,1:biter]

for (k in 1:kappa.cont.n) {
  boot.cont.mean <- apply(pred.cont.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.cont.sd <- apply(pred.cont.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.cont.update[,,z] <- ifelse((pred.cont.update[,,z] < (boot.cont.mean-kappa*boot.cont.sd) | pred.cont.update[,,z] >
                                  (boot.cont.mean+kappa*boot.cont.sd)), NA, pred.cont.update[,,z])
  } # end z
  print(k)
} # end k

pred.cont.med <- apply(pred.cont.update, MARGIN=c(1,2), median, na.rm=T)
pred.cont.lo <- apply(pred.cont.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.cont.up <- apply(pred.cont.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.cont.med <- apply(cont.val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

## factor variables
# gender
pred.fac1.update <- fac1.pred.mat[1:biter,]
boot.fac1.mean <- apply(fac1.pred.mat, MARGIN=2, mean, na.rm=T)
boot.fac1.sd <- apply(fac1.pred.mat, MARGIN=2, sd, na.rm=T)

for (z in 1:biter) {
  pred.fac1.update[z, ] <- ifelse((fac1.pred.mat[z, ] < (boot.fac1.mean-kappa*boot.fac1.sd) | fac1.pred.mat[z, ] >
                                  (boot.fac1.mean+kappa*boot.fac1.sd)), NA, fac1.pred.mat[z, ])
} # end z

# frailty
boot.fac2.mean <- apply(fac2.pred.mat, MARGIN=2, mean, na.rm=T)
boot.fac2.sd <- apply(fac2.pred.mat, MARGIN=2, sd, na.rm=T)

for (z in 1:biter) {
  fac2.pred.mat[z, ] <- ifelse((fac2.pred.mat[z, ] < (boot.fac2.mean-kappa*boot.fac2.sd) | fac2.pred.mat[z, ] >
                                  (boot.fac2.mean+kappa*boot.fac2.sd)), NA, fac2.pred.mat[z, ])
} # end z

pred.fac1.med <- apply(fac1.pred.mat, MARGIN=2, median, na.rm=T)
pred.fac1.lo <- apply(fac1.pred.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
pred.fac1.up <- apply(fac1.pred.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
pred.fac2.med <- apply(fac2.pred.mat, MARGIN=2, median, na.rm=T)
pred.fac2.lo <- apply(fac2.pred.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
pred.fac2.up <- apply(fac2.pred.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:npredictcols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(predictor_vars,".mean",sep="")
vec.sd.names <- paste(predictor_vars,".sd",sep="")

for (k in 1:npredictcols) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:npredictcols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:npredictcols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med <- median(CV.cor.update, na.rm=T)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo, CV.cor.med, CV.cor.up))

ri.vec.lo.names <- paste(predictor_cols,".ri.lo",sep="")
ri.vec.up.names <- paste(predictor_cols,".ri.up",sep="")
ri.vec.med.names <- paste(predictor_cols,".ri.med",sep="")

for (ri in 1:npredictcols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- predictor_vars
ri.sort <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort

# plot
ri.plt <- ggplot(ri.sort) +
  geom_bar(aes(x=reorder(row.names(ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo,CV.cor.med,CV.cor.up), 2))

## plot predicted relationships
# continuous variables
head(pred.cont.med)
ri.sort
pred.cont.bs.med <- logit2prob(pred.cont.med)
pred.cont.bs.lo <- logit2prob(pred.cont.lo)
pred.cont.bs.up <- logit2prob(pred.cont.up)
age.orig <- RESP.cont.val.dat$age.sc * attr(infect.dat2$age.sc, "scaled:scale") + attr(infect.dat2$age.sc, "scaled:center")
bmi.orig <- 10^(RESP.cont.val.dat$lbmi.sc * attr(infect.dat2$lbmi.sc, "scaled:scale") + attr(infect.dat2$lbmi.sc, "scaled:center"))
predictors.orig <- as.data.frame(cbind(age.orig, bmi.orig))
xlabs.plot <- c("patient age (years)", expression("patient body mass index (kg/m"^2*")"))
ylims <- c(min(pred.cont.bs.lo[,cont_predictors], na.rm=T), max(pred.cont.bs.up[,cont_predictors], na.rm=T))
plotNcontvec <- paste("plt",1:length(cont_predictor_cols),sep="")
for (v in 1:length(cont_predictor_cols)) {
  dat.use <- as.data.frame(cbind(predictors.orig[v], pred.cont.bs.med[,cont_predictors[v]],
                                 pred.cont.bs.lo[,cont_predictors[v]], pred.cont.bs.up[,cont_predictors[v]]))
  colnames(dat.use) <- c("x", "V2", "V3", "V4")
  assign(plotNcontvec[v], ggplot(data=dat.use) +
           geom_line(aes(x=x, y=V2), colour="blue") +
           geom_ribbon(aes(x=x, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(xlabs.plot[v]) + ylab("Pr(infection)"))
}
ggarrange(plt1, plt2, ncol=2, nrow=1)

# factor variables (gender)
pred.bs.fac1.med <- logit2prob(pred.fac1.med)
pred.bs.fac1.lo <- logit2prob(pred.fac1.lo)
pred.bs.fac1.up <- logit2prob(pred.fac1.up)
plt3.dat <- as.data.frame(cbind(gender=rownames(RESP.fac1.pred.dat), med=pred.bs.fac1.med, lo=pred.bs.fac1.lo, up=pred.bs.fac1.up))
plt3.dat$med <- format(signif(as.numeric(plt3.dat$med), digits = 3), trim = TRUE)
plt3.dat$lo <- format(signif(as.numeric(plt3.dat$lo), digits = 3), trim = TRUE)
plt3.dat$up <- format(signif(as.numeric(plt3.dat$up), digits = 3), trim = TRUE)

plt3 <- ggplot(data=plt3.dat) +
  geom_bar(aes(x=gender, y=med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=gender, ymin=lo, ymax=up),
                linewidth=0.4, colour="black", alpha=0.9) +
  xlab("patient gender") + ylab("Pr(infection)")

# factor variables (frailty)
pred.bs.fac2.med <- logit2prob(pred.fac2.med)
pred.bs.fac2.lo <- logit2prob(pred.fac2.lo)
pred.bs.fac2.up <- logit2prob(pred.fac2.up)
plt4.dat <- as.data.frame(cbind(frailty=rownames(RESP.fac2.pred.dat), med=pred.bs.fac2.med, lo=pred.bs.fac2.lo, up=pred.bs.fac2.up))
plt4.dat$med <- format(signif(as.numeric(plt4.dat$med), digits = 3), trim = TRUE)
plt4.dat$lo <- format(signif(as.numeric(plt4.dat$lo), digits = 3), trim = TRUE)
plt4.dat$up <- format(signif(as.numeric(plt4.dat$up), digits = 3), trim = TRUE)

plt4 <- ggplot(data=plt4.dat) +
  geom_bar(aes(x=frailty, y=med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=frailty, ymin=lo, ymax=up),
                linewidth=0.4, colour="black", alpha=0.9) +
  xlab("patient frailty (1=fit ... 8=extremely frail)") + ylab("Pr(infection)")

ggarrange(plt1, plt2, plt3, plt4, ncol=2, nrow=2)

