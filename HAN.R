


getwd()
setwd("/Users/uisupshin/Dropbox/HANJM")
getwd()

# Call the package
#install moonBook from github
#install.packages("devtools")  # devtools
#install.packages("moonBook")
require(moonBook)
require(survival)

#install.packages("maxstat")
require(maxstat)
require(grid)
#install.packages("ggplot2")
require(ggplot2)
source("ggKM.R")
source("mycphwt.R")

HAN_raw=read.csv("han20230213.csv", header=T, dec=".")

# Select Rows by equal condition
HAN <- subset(HAN_raw, Adjuvant_Doublet==1)

# Variable Set up
HAN$LNR<-HAN$pN/HAN$tN
HAN$PD<-ifelse(HAN$Diff=="PD",1,0)
HAN$Pni1<-ifelse(HAN$N=="Y",1,0)
HAN$Vi1<-ifelse(HAN$V=="Y",1,0)
HAN$Lvi1<-ifelse(HAN$L=="Y",1,0)
HAN$Male<-ifelse(HAN$SEX==1,1,0)
HAN$RTCOL<-ifelse(HAN$SIDENESS=="Rt",1,0)
HAN$CCI7<-ifelse(HAN$CCI>7,1,0)
HAN$AGE60<-ifelse(HAN$AGE>60,1,0)
HAN$CEA73<-ifelse(HAN$preHR_CEA>73,1,0)
HAN$ASA3<-ifelse(HAN$ASA>2,1,0)
HAN$Adju_target<-ifelse(HAN$Post_HR_target_all=="None",0,ifelse(HAN$Post_HR_target_all=="Avastin",1,2))
HAN$Adju_target <-as.factor(HAN$Adju_target)

# Time object, DFS
HAN$TS_DFS=Surv(HAN$FU_DFS, HAN$RECUR=="1") #time object
fit_DFS <-survfit(HAN$TS_DFS~., data=HAN)



#Cut-off for DFS
####CCI
maxstat.test(TS_DFS~ preHR_CEA, data=HAN,smethod="LogRank",pmethod="condMC",B=999)
maxstat.test(TS_DFS~ CCI, data=HAN,smethod="LogRank",pmethod="condMC",B=999)
maxstat.test(TS_DFS~ ASA, data=HAN,smethod="LogRank",pmethod="condMC",B=999)
maxstat.test(TS_DFS~ NUM_HM, data=HAN,smethod="LogRank",pmethod="condMC",B=999)
maxstat.test(TS_DFS~ HM_SIZE, data=HAN,smethod="LogRank",pmethod="condMC",B=999)
maxstat.test(TS_DFS~ LNR, data=HAN,smethod="LogRank",pmethod="condMC",B=999)

HAN$CCI7<-ifelse(HAN$CCI>7,"CCI>7","CCI<7")
DFS<-survfit(Surv(FU_DFS, RECUR=="1")~Adju_TA, data=HAN)
summary(DFS)
DFS

plot(DFS, main="Progression free survival, Target agent, yes/no", xlab="Time", ylab="Survival Probability", col=c("blue", "red"), lty=c(1,2))
legend("topright", legend=c("TA no","TA yes"), col=c("blue", "red"), lty=c(1,2), cex=0.8)

head(HAN)


table1<-mytable(Adju_TA~AGE+Male+Synch_M1+M1_12m+ASA3+CCI7+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+preHR_CEA+CEA73+PreHR_chemo+PreHR_TA+Post_HR_target_all,
                show.total=TRUE, exact=TRUE, data=HAN, method = 3)
table1

table2<-mytable(Adju_target~AGE+Male+Synch_M1+M1_12m+ASA3+CCI7+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+preHR_CEA+CEA73+PreHR_chemo+PreHR_TA+Post_HR_target_all,
                show.total=TRUE, exact=TRUE, data=HAN, method = 3)
table2

mycsv(table1, file="table1.csv")


# Time object, DFS
HAN$TS_DFS=Surv(HAN$FU_DFS, HAN$RECUR=="1") #time object
fit_DFS <-survfit(HAN$TS_DFS~Adju_TA, data=HAN)

### UNIVARIATE
out_u=mycph(TS_DFS~Adju_TA+Any_TA+Adju_target+Male+AGE60+CCI7+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+CEA73+M1_12m+PreHR_chemo+PreHR_TA, data=HAN) # Univariate
out_u
write.csv(out_u, file="Cox_uni_DFS.csv")

#MULTI, p=0.25
out_m<-coxph(TS_DFS~Any_TA+Male+AGE60+Vi1+BILOBE+MULTI4+CEA73, data=HAN)

extractHR(out_m)
DFS_step<-step(out_m, direction="backward")
Cox_step_dfs<-extractHR(DFS_step)
Cox_step_dfs
write.csv(Cox_step_dfs, file="Cox_step_DFS.csv")

HAN$TS_OS=Surv(HAN$FU_OS, HAN$DEATH==1) #time object
fit_OS <-survfit(HAN$TS_OS~Adju_TA, data=HAN)

out_u_os=mycph(TS_OS~Adju_TA+SEX+AGE60+CCI7+RTCOL++LNR+pT4+Pni1+Lvi1+Vi1+PD, data=HAN) # Univariate
out_u_os

#MULTI_os
out_m_os<-coxph(TS_OS~Adju_TA+SEX+AGE60+CCI7+RTCOL++LNR+pT4+Pni1+Lvi1+Vi1+PD, data=HAN)

extractHR(out_m_os)




## Survival Curve
#install.packages("survminer")
library(survminer)
DFS_all<-survfit(Surv(FU_DFS, RECUR=="1")~1, data=HAN)

summary(DFS_all)
min(HAN$FU_DFS)
max(HAN$FU_DFS)
median(HAN$FU_DFS)
median(HAN$FU_OS)
min(HAN$FU_OS)
max(HAN$FU_OS)
sum(HAN$RECUR)
sum(HAN$DEATH)

DFS<-survfit(Surv(FU_DFS, RECUR=="1")~Adju_TA, data=HAN)
summary(DFS)
DFS

ç<-survfit(Surv(FU_DFS, RECUR=="1")~Adju_target, data=HAN)
summary(DFS_target_kind)

plot(DFS_target_kind, main="Progression free survival, Kind of target agent", xlab="Time", ylab="Survival Probability", col=c("blue", "red","green"), lty=c(1,2))
legend("topright", legend=c("None","Bevacizumab","Cetuximab"), col=c("blue", "red","green"), lty=c(1,2), cex=0.8)

ggsurvplot(DFS_target_kind, 
           legend=c(0.8,0.8),
           pval=TRUE, 
           xlim=c(0,60),
           break.time.by = 12,
           legend.labs = 
             c("None", "Bevacizumab","Cetuximab"))


ggsurvplot(
  DFS,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,        # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 12,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  #surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Adjuvant target agent, No", "Yes"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF"), # custom color palettes.
  fun = NULL, #Cumulative risk
  xlim = c(0, 72)
)


#### Overall Survival ####
OS<-survfit(Surv(FU_OS, DEATH=="1")~Adju_TA, data=HAN)
OS_all<-survfit(Surv(FU_OS, DEATH=="1")~1, data=HAN)
summary(OS_all)
OS

ggsurvplot(
  OS,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,        # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 12,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  #surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Adjuvant target agent, No", "Yes"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF"), # custom color palettes.
  fun = NULL, #Cumulative risk
  xlim = c(0, 72)
)
############### IPTW ####################################################

library(dplyr)
library(MatchIt)

HAN_IPTW<-as.data.frame(HAN[,c("MRN","Adju_TA","Adju_target","Male","AGE60","CCI7","ASA3","RTCOL","LNR","pT4","Pni1","Lvi1","Vi1","PD","BILOBE","MULTI4","SIZE4"
                               ,"CEA73","M1_12m","PreHR_chemo","PreHR_TA","FU_DFS","RECUR")])


## DATA character change ###
HAN_IPTW$Male <-as.factor(HAN_IPTW$Male)
HAN_IPTW$AGE60 <-as.factor(HAN_IPTW$AGE60)
HAN_IPTW$CCI7 <-as.factor(HAN_IPTW$CCI7)
HAN_IPTW$ASA3 <-as.factor(HAN_IPTW$CCI7)
HAN_IPTW$RTCOL <-as.factor(HAN_IPTW$RTCOL)
HAN_IPTW$pT4 <-as.factor(HAN_IPTW$pT4)
HAN_IPTW$Pni1 <-as.factor(HAN_IPTW$Pni1)
HAN_IPTW$Lvi1 <-as.factor(HAN_IPTW$Lvi1)
HAN_IPTW$Vi1 <-as.factor(HAN_IPTW$Vi1)
HAN_IPTW$PD <-as.factor(HAN_IPTW$PD)
HAN_IPTW$BILOBE <-as.factor(HAN_IPTW$BILOBE)
HAN_IPTW$MULTI4 <-as.factor(HAN_IPTW$MULTI4)
HAN_IPTW$SIZE4 <-as.factor(HAN_IPTW$SIZE4)
HAN_IPTW$CEA73 <-as.factor(HAN_IPTW$CEA73)
HAN_IPTW$M1_12m <-as.factor(HAN_IPTW$M1_12m)
HAN_IPTW$PreHR_TA <-as.factor(HAN_IPTW$PreHR_TA)



table1_iptw <- mytable(Adju_TA~Adju_target+Male+AGE60+CCI7+ASA3+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+CEA73+M1_12m+PreHR_TA,
                       show.total=TRUE, exact=TRUE, data=HAN_IPTW, method = 3)

table1_iptw

## Calculating propensity score ###

fit_IPTW<-glm(Adju_TA~Male+AGE60+ASA3+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+CEA73+M1_12m+PreHR_TA, family=binomial, data=HAN_IPTW)
factor<-extractOR(fit_IPTW, digits=3)
factor
summary(fit_IPTW)

write.csv(factor_PTR, file="fit_PTR.csv")

HAN_IPTW$pr.score <- predict(fit_IPTW, type = "response") ##attach propensity score##
#HAN_IPTW$att.weights <- with(HAN_IPTW, Adju_TA + (1-Adju_TA)*pr.score/(1-pr.score))
HAN_IPTW$ate.weights <- ifelse(HAN_IPTW$Adju_TA=="1", 1/HAN_IPTW$pr.score, 1/(1-HAN_IPTW$pr.score))


## Goodness of Fit test, Hosmer Lemeshow test ####
#install.packages("ResourceSelection")
library(ResourceSelection)
hoslem.test(HAN_IPTW$Adju_TA, fitted(fit_IPTW))

library(pROC)
#install.packages("Epi")
require(Epi)
#install.packages("pROC")
a1=ROC(form=Adju_TA~pr.score, data=HAN_IPTW, plot="ROC") ## discrimination of TA NO vs. TA YES group by propensity score ##

## Weightit packages ###
#install.packages("WeightIt")
library (WeightIt)
library(cobalt)
covs <- subset(HAN_IPTW, select = c(Male,AGE60,ASA3,RTCOL,LNR,pT4,Pni1,Lvi1,Vi1,PD,BILOBE,MULTI4,SIZE4,CEA73,M1_12m,PreHR_TA))

#W.out_att <- weightit(Adju_TA~covs, data=HAN_IPTW, method="ps", estimand = "ATT") ##table#
#bal.tab(W.out_att, un=TRUE, binary="std", thresholds=0.2)

W.out_ate <- weightit(Adju_TA~covs, data=HAN_IPTW, method="ps", estimand = "ATE") ##table#
bal.tab(W.out_ate, un=TRUE, binary="std", thresholds=0.2)

## Love plot 1

love.plot(W.out_ate, thresholds = c(m = .1), var.order = "unadjusted")


## Love plot 2, # balance plot#, Figure 1##
v<- data.frame(old=c("Male","AGE60","ASA3","RTCOL","LNR","pT4","Pni1","Lvi1","Vi1","PD","BILOBE","MULTI4","SIZE4","CEA73","M1_12m","PreHR_TA","ps"),
               new=c("Male","Age>60","ASA >3","Primary tumor at Rt colon","Regional LN ratio","pT4", "Perineural invasion of PT",
                    "Lymphatic invasion of PT","Vascular invasion of PT", "Poor differentiated", "Bilobar involvement of CRLM", 
                    "Number of CRLM>4", "Max size of CRLM>4cm", "Pre-hepatectomy CEA>73", "Synch or LM within 12mon","Neoaduvant target agent", "propensity score"))

love.plot(W.out_ate, stats = c("mean.diffs"), 
          threshold = c(m = .2), 
          binary = "std", abs = TRUE,
          var.order = "unadjusted",
          var.names = v,
          limits = c(0, 1), grid = TRUE, wrap = 20,
          sample.names = c("Original", "Weighted"),
          position = "top",
          line= TRUE,
          shapes = c("circle", "triangle"),
          colors = c("red", "blue"))+
  labs(x="Absolute standadized mean difference")+
  theme(legend.position = c(.75, .3),
        legend.box.background = element_rect(), 
        legend.box.margin = margin(0.5, 1, 1, 1),
        legend.title=element_blank(),
        plot.title=element_blank(),
        axis.title=element_text(color = "black", size=10))

### TableOne package ####
#install.packages("tableone")
library(tableone)
#install.packages("survey")
library(survey)

vars <- c("Male","AGE60","ASA3","RTCOL","LNR","pT4","Pni1","Lvi1","Vi1","PD","BILOBE","MULTI4","SIZE4","CEA73","M1_12m","PreHR_TA")

tab.unwt <- CreateTableOne(vars=vars, strata="Adju_TA", data=HAN_IPTW, test =FALSE)
print(tab.unwt, smd=TRUE)

Svy <- svydesign(ids = ~ 1, data = HAN_IPTW, weights = ~ ate.weights)
summary(Svy)
Svy
tab.wt <- svyCreateTableOne(vars=vars, strata="Adju_TA", data=Svy, test =FALSE)
print(tab.wt, smd=TRUE)



library(jskm)   

unadjusteDFS<-survfit(Surv(HAN_IPTW$FU_DFS, as.numeric(HAN_IPTW$RECUR)) ~ HAN_IPTW$Adju_TA, data=HAN_IPTW)
summary(unadjusteDFS)

iptwsvyDFS <- svydesign(ids = ~ 1, data = HAN_IPTW, weights = ~ ate.weights)

iptwDFS<-svykm(Surv(FU_DFS,RECUR)~Adju_TA, design=iptwsvyDFS, data=HAN_IPTW)

iptwDFS

a<-jskm(unadjusteDFS, timeby = 12, ystratalabs=c('No',"Yes"),
        ystrataname = "Adjuvant Target therapy", table = TRUE, ci=FALSE, pval=FALSE,
        xlabs="Months after Surgery", main ="Unadjusted",
        dashed=FALSE,marks=FALSE,xlims=c(0,60),legendposition=c(0.85,0.85))

a


b<-svyjskm(iptwDFS, timeby = 12, ystratalabs=c('No',"Yes"),
           ystrataname = "Adjuvant Target therapy", table = TRUE, ci=FALSE, pval=FALSE,
           xlabs="Months after Surgery", main ="Weighted",
           dashed=FALSE,marks=TRUE,xlims=c(0,60),legendposition=c(0.85,0.85), pval=TRUE)

b

#### Univarate Cox regression, unweighted, weighted ######################

HAN_IPTW$TS_DFS=Surv(HAN_IPTW$FU_DFS, HAN_IPTW$RECUR=="1") #time object

## Uni, unweighted
out = mycph(TS_DFS~Adju_TA+Adju_target+Male+AGE60+ASA3+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+CEA73+M1_12m+PreHR_TA, digits=2, data=HAN_IPTW)
out
write.csv(out, file="Cox_uni.csv")


## Uni, weighted
uni_ate = mycphwt(TS_DFS~Adju_TA+Adju_target+Male+AGE60+ASA3+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+CEA73+M1_12m+PreHR_TA, digits=2, data=HAN_IPTW)
uni_ate
write.csv(uni_ate, file="Cox_uni_ate.csv")

##### Multivariate Cox regression, unweighted, weighted #################

## unweighted multi cox
COX<-coxph(TS_DFS~Adju_TA+Vi1+BILOBE+MULTI4+CEA73+PreHR_TA, data=HAN_IPTW)
COX<-coxph(TS_DFS~Adju_target+Vi1+BILOBE+MULTI4+CEA73+PreHR_TA, data=HAN_IPTW)
summary(COX)
extractHR(COX)

COX_step<-step(COX, direction="backward")
extractHR(COX_step, digit=2)

## weighted multi cox###
## ATE ###
multi_ate<-coxph(TS_DFS~Adju_TA+Male+LNR+Pni1+Vi1+BILOBE+MULTI4+CEA73, weights=ate.weights, data=HAN_IPTW)
multi_ate<-coxph(TS_DFS~Adju_target+Male+LNR+Pni1+Vi1+BILOBE+MULTI4+CEA73, weights=ate.weights, data=HAN_IPTW)
summary(multi_ate)

extractHRwt(multi_ate)

Wt_step<-step(multi_ate, direction="backward")
summary(Wt_step)
t.step<-extractHRwt(Wt_step, digit=2)
t.step
write.csv(t.step, file="multi_Cox_wt.csv")


##Any_TA+Male+AGE60+CCI7+RTCOL+LNR+pT4+Pni1+Lvi1+Vi1+PD+BILOBE+MULTI4+SIZE4+CEA73+M1_12m+PreHR_chemo+PreHR_TA