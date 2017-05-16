#########################################################################################
################## This file imports the data and cleans it for analysis  ###############
#########################################################################################

# the cleaned data set used for analysis is stored as lymph.rda
lymph <- read.delim(file = '~/Repositories/Data/Causal_Analytics/Lymphom.txt', sep = '')

# have a look
library(psych)
describe(lymph)
# no missing, imported correctly

# relabel sex 
levels(lymph$Sex) <- list('Female' = 'F', 'Male' = 'M')

# factor categorical variables
lymph$Reg <- factor(lymph$Reg)
lymph$Rad <- factor(lymph$Rad, labels = c('No', 'Yes'))
lymph$RdSrg <- factor(lymph$RdSrg, labels = c("No", "Yes"))
lymph$Site <- factor(lymph$Site, labels = c('non-Hodgkin', 'Hodgkin'))
lymph$RacB <- factor(lymph$RacB, labels = c('White', 'Black'))
lymph$Stag <- factor(lymph$Stag, labels = c('Localized', 'Regional', 
                                            'Distant', 'Unknown'))
lymph$Caus <- factor(lymph$Caus, labels = c('Alive', 'Death From Lymphoma', 
                                            'Death from non-lymphoma Cancer', 
                                            'Death from non-cancer cause'))
lymph$LymDth <- lymph$Caus == 'Death From Lymphoma'

# look at death and cause of death
table(lymph$Dth, lymph$Caus, useNA = 'ifany')
# They match!
table(lymph$Caus, lymph$LymDth, useNA = 'ifany')
# great

# will exclude any one who did not recieve radiation as they are not relevant to our analyis
# will store number excluded for final report
# check to make sure there are no 'no radiation' plus 'yes to radiation and surgery'
table(lymph$Rad, lymph$RdSrg, useNA = 'ifany')
# checks out!
no_rad <- nrow(lymph[lymph$Rad == 'No',])

#remove no radiation
lymph <- lymph[lymph$Rad == 'Yes',]
# 22505 removed, 9184 remaining

save(lymph, file = '~/Repositories/Data/Causal_Analytics/lymph.rda')
#########################################################################################
########################### look at descriptives/distributions  #########################
#########################################################################################
boxplot(Age ~ RdSrg, data = lymph)
boxplot(Grad ~ RdSrg, data = lymph)
plot(Reg ~ RdSrg, data = lymph)
plot(Dth ~ RdSrg, data = lymph)
plot(Sex ~ RdSrg, data = lymph)
plot(Site ~ RdSrg, data = lymph)
plot(RacB ~ RdSrg, data = lymph)
plot(Stag ~ RdSrg, data = lymph)
boxplot(Tim ~ RdSrg, data = lymph)

names(lymph)[c(1:2,5,10:11)] <- c('Registry/Location', 'Birth Year', 'Grade of Lymphoma', 'Race', 
                                 'SEER Stage')
source('~/Repositories/Table1/Table1.R')
tbl1unw <- Table1(c(1:5,9:11), 7, lymph)




save.image(file = '~/Repositories/Data/Causal_Analytics/dataclean.rdata')
