################################################################################
######### This code contains all code used in the statisitcal analysis #########
################################################################################

# for all analysis I used guidelines from: 
# Austin PC. The use of propensity score methods with survival or 
# time-to-event outcome: reporting measures of effect similar to 
# those used in randomized experiments. Statistics in Medicine 2014;
# 33:1242-1258

library(survey)
library(survival)

load(file = '~/Repositories/Data/Causal_Analytics/lymph.rda')

################################################################################
############################### Propensity Score ###############################
propscore <- glm(RdSrg ~ Reg + BYr + Age + Sex + Grad + Site + 
                   RacB + Stag, data = lymph, family = binomial())
fitted.probs <- predict(propscore, type = 'response')
lymph$PropScore <- fitted.probs


# Inverse Probability of Treatment Weighting Using the Propesity Score
w1=ifelse(lymph$RdSrg=="Yes", 1/fitted.probs, 1/(1 - fitted.probs))
design=svydesign(ids=~1,weights=~w1,data=lymph)

# analyze balance
#get variable names from dataser
rv <- names(lymph)[c(1:5,9:11)]
source(file = '~/Repositories/Table1/Table1Weighted.R')
Tbl1wgt <- Table1Weighted(rv, 'RdSrg', design)
#nice!


################################################################################
############################# kaplan meier curves ##############################
# crude 
km_crude <- survfit(Surv(Tim, Dth) ~ RdSrg, data = lymph, conf.type = "log-log")
plot(km_crude, lwd = 2, lty = c(1,2), col = c("purple","forestgreen"), 
     xlab = 'Time to Death (days)', 
     ylab = 'Survival Probability')
legend(120, 1, c('Radiation w/o Surgery', 'Radiation with Surgery'), lty = 1:2, 
       col=c("purple","forestgreen"), bty = "n")
title("Crude Kaplan-Meier Curve")

# weighted
km_wgt <- svykm(Surv(Tim, Dth) ~ RdSrg, design = design)
km_wgt

plot(km_wgt, lwd=2, pars=list(lty=c(1,2),col=c("purple","forestgreen")), 
     xlab = 'Time to Death (days)', 
     ylab = 'Survival Probability')
legend(120, 1, c('Radiation w/o Surgery', 'Radiation with Surgery'), lty = 1:2, 
       col=c("purple","forestgreen"), bty = "n")
title('Inverse Probability Weighted KM Curve')

################################################################################
############################# cox hazard regression ############################

model <- svycoxph(Surv(Tim, Dth) ~ RdSrg, design = design)
model

model_lymdth <- svycoxph(Surv(Tim, LymDth) ~ RdSrg, design = design)
model_lymdth

################################################################################
############################# Sensitivity Analysis #############################
library(obsSens)

sensitivity <- obsSensSCC(model, which = 1, 
                          g0 = c(2, 6, 10), 
                          p0 = seq(0.01,0.4,by=0.05), 
                          p1 = seq(0.01,0.4,by=0.05), 
                          logHaz = FALSE, method = c("approx", "sim"))
summary(sensitivity)

save.image(file = '~/Repositories/Data/Causal_Analytics/analysis.rdata')
