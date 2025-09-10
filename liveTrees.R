library(dplyr)
library(ggplot2)


treeData <- read.csv("tree_data.csv", header = T)

## Data from 2013, example reference date is 1960 (53 years ago)

treeData <- treeData %>% 
  select(Transect..:DBH)

treeData$DBH <- as.numeric(treeData$DBH)
treeData$Age <- as.numeric(treeData$Age)


lm <- lm(data = treeData, formula = Age ~ DBH)

summary(lm)

coefficients(lm)
coefficients(lm)[1]
coefficients(lm)[2]

AgeDBHplot <- ggplot(data = treeData, mapping = aes(x = DBH, y = Age)) +
  geom_point() +
  geom_smooth(method = "lm")

AgeDBHplot


predictedAge <- predict(lm, newdata = treeData)

ageTreeData <- treeData

ageTreeData$Age[is.na(ageTreeData$Age)] <- predictedAge[is.na(ageTreeData$Age)]
ageTreeData$Age <- ceiling(ageTreeData$Age)


ageTreeData$estabYear <- 2013 - ageTreeData$Age


refTreeData <- ageTreeData %>% 
  filter(estabYear <= 1960) %>% 
  mutate(DBHrefYear = DBH - ((2013-1960)*(coefficients(lm)[2]))) # stuck right here




