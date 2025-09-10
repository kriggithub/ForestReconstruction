library(dplyr)
library(ggplot2)


treeData <- read.csv("tree_data.csv", header = T)

## Data from 2013, example reference date is 1960 (53 years ago)

treeData <- treeData %>% 
  select(Transect..:DBH)

treeData$DBH <- as.numeric(treeData$DBH)
treeData$Age <- as.numeric(treeData$Age)
treeData$Species <- trimws(treeData$Species)

# or use treeData$Species[treeData$Species == "ABBI "] <- "ABBI"

treeData$Species <- toupper(treeData$Species)


unique(treeData$Species)


# Code not for by species
# lm <- lm(data = treeData, formula = Age ~ DBH)
# 
# summary(lm)
# coefficients(lm)
# coefficients(lm)[1]
# coefficients(lm)[2]
# 
# 
# AgeDBHplot <- ggplot(data = treeData, mapping = aes(x = DBH, y = Age)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# AgeDBHplot




######################
# For trees that are alive


# linear models per species to determine non-aged trees

lmABBI <- lm(data = treeData, formula = Age ~ DBH, subset = Species == "ABBI")
lmPIEN <- lm(data = treeData, formula = Age ~ DBH, subset = Species == "PIEN")


summary(lmABBI)
coefficients(lmABBI)
coefficients(lmABBI)[1]
coefficients(lmABBI)[2]

summary(lmPIEN)
coefficients(lmPIEN)
coefficients(lmPIEN)[1]
coefficients(lmPIEN)[2]



AgeDBHplotABBI <- ggplot(data = treeData, mapping = aes(x = DBH, y = Age)) +
  geom_point(data = subset(treeData, Species == "ABBI")) +
  geom_smooth(data = subset(treeData, Species == "ABBI"), method = "lm")

AgeDBHplotABBI


AgeDBHplotPIEN <- ggplot(data = treeData, mapping = aes(x = DBH, y = Age)) +
  geom_point(data = subset(treeData, Species == "PIEN")) +
  geom_smooth(data = subset(treeData, Species == "PIEN"), method = "lm")

AgeDBHplotPIEN


# Not species specific
# predictedAge <- predict(lm, newdata = treeData)
# 
# ageTreeData <- treeData
# 
# ageTreeData$Age[is.na(ageTreeData$Age)] <- predictedAge[is.na(ageTreeData$Age)]
# ageTreeData$Age <- ceiling(ageTreeData$Age)




# create new dataframe with new rows, filling depending on species

ageTreeData <- treeData
missingRows <- is.na(ageTreeData$Age)
ABBIrows <- missingRows & ageTreeData$Species == "ABBI"


if (any(ABBIrows)) {
  predictedAgeABBI <- predict(lmABBI, newdata = ageTreeData[ABBIrows, ])
  ageTreeData$Age[ABBIrows] <- predictedAgeABBI
}



PIENrows <- missingRows & ageTreeData$Species == "PIEN"


if (any(PIENrows)) {
  predictedAgePIEN <- predict(lmPIEN, newdata = ageTreeData[PIENrows, ])
  ageTreeData$Age[PIENrows] <- predictedAgePIEN
}




# round age up to next whole number
ageTreeData$Age <- ceiling(ageTreeData$Age)



ageTreeData$estabYear <- 2013 - ageTreeData$Age


###################################
######### NOTE
###################################
# Datasheets for RMBL Reconstruction suggest that Avg incr. (mm) seem to be per transect
# Question: How was this determined?

# For now, assume ABBI avg. increment = 0.822818 mm/year
# PIEN avg. increemnt = 0.85 mm/year

ABBIavgInc <- 0.822818
PIENavgInc <- 0.85

# Then, growth to subtract for ABBI & PIEN (cm per year) = 

ABBIgrowthSub <- (ABBIavgInc*2)/10
PIENgrowthSub <- (PIENavgInc*2)/10



refTreeData <- ageTreeData %>% 
  filter(estabYear <= 1960) %>% 
  mutate(RefAge = (1960 - estabYear),
         RefDBH = ifelse(
           Species == "ABBI",
           DBH - ((2013 - 1960) * ABBIavgInc),
           ifelse(
             Species == "PIEN",
             DBH - ((2013 - 1960) * PIENavgInc),
             NA))) %>% 
  filter(RefDBH >= 0)


# Live tree DBHs calculated






############################# 
# Reviving dead trees to their state back at the reference date


# Calculate avg. dbh for species and convert to inches


# Don't know where this came from but avg. in dataset is (cm):
ABBIavgDBH <- 14.8
PIENavgDBH <- 19.1
ABBIavgDBHin <- 14.8 * 2.54
PIENavgDBHin <- 19.1 * 2.54



# Snag Life  for (GUMO) = 2 * DBH(in inches)
# Snag Fall = 1/Snag Life

ABBIsl <- 2 * ABBIavgDBHin
PIENsl <- 2 * PIENavgDBHin

ABBIsf <- 1/ABBIsl
PIENsf <- 1/PIENsl




# decoposition rates between conclasses
# 3 -> 4 = 20%/year
# 4 -> 5 = 20%/year
# 5 -> 6 = 20%/year
# 6 -> 7 = 20%/year
# 7 -> 8 = 20%/year




























save.image(file = "liveTrres.RData")

