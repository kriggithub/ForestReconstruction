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
ABBIavgDBHin <- ABBIavgDBH/2.54
PIENavgDBHin <- PIENavgDBH/2.54



# Snag Life  for (GUMO) = 2 * DBH(in inches)
# Snag Fall = 1/Snag Life

ABBIsl <- 2 * ABBIavgDBHin
PIENsl <- 2 * PIENavgDBHin

ABBIsf <- 1/ABBIsl
PIENsf <- 1/PIENsl




# decomposition rates between conclasses

speciesnames <- c("ABBI", "PIEN")
conclasses <- 3:7


decompreference <- expand.grid(Species = speciesnames,
                               Conclass = conclasses)

decompreference$Rate <- NA

decompreference$Rate[decompreference$Species == "ABBI" & decompreference$Conclass == 3] <- 0
decompreference$Rate[decompreference$Species == "PIEN" & decompreference$Conclass == 3] <- 0
decompreference$Rate[decompreference$Species == "ABBI" & decompreference$Conclass == 4] <- 0.2
decompreference$Rate[decompreference$Species == "PIEN" & decompreference$Conclass == 4] <- 0.2
decompreference$Rate[decompreference$Species == "ABBI" & decompreference$Conclass == 5] <- 0.15
decompreference$Rate[decompreference$Species == "PIEN" & decompreference$Conclass == 5] <- 0.15
decompreference$Rate[decompreference$Species == "ABBI" & decompreference$Conclass == 6] <- ABBIsf
decompreference$Rate[decompreference$Species == "PIEN" & decompreference$Conclass == 6] <- PIENsf
decompreference$Rate[decompreference$Species == "ABBI" & decompreference$Conclass == 7] <- ABBIsf
decompreference$Rate[decompreference$Species == "PIEN" & decompreference$Conclass == 7] <- PIENsf
# decompreference$Rate[decompreference$Species == "ABBI" & decompreference$Conclass == "8"] <- ABBIsf
# decompreference$Rate[decompreference$Species == "PIEN" & decompreference$Conclass == "8"] <- PIENsf






percentiles <- c(0.25, 0.5, 0.75)

decompRate <- function(rate, percentiles){
    (log(percentiles) - log(1))/(log(1 + rate))
}


for (p in percentiles){
  colname <- as.character(p)
  
  unaddedDecompRate <- decompRate(decompreference$Rate, p)
  
  decompreference[[colname]] <- ave(unaddedDecompRate, decompreference$Species, FUN = cumsum)

}









percentiles <- c(0.25, 0.5, 0.75)

# 3 -> 4 = 20%/year (known)
rate34ALL <- 0.2
decomprate34ALL <- decompRate(rate34ALL, percentiles)


# 4 -> 5 = 15%/year(known)
rate45ALL <- 0.15
decomprate45ALL <- decompRate(rate45ALL, percentiles)
decomprate45ALL <- decomprate34ALL + decomprate45ALL

# 5 -> 6 = dependent on snag fall
rate56PIEN <- PIENsf
decomprate56PIEN <- decompRate(rate56PIEN,  percentiles)
decomprate56PIEN <- decomprate56PIEN + decomprate45ALL

rate56ABBI <- ABBIsf
decomprate56ABBI <- decompRate(rate56ABBI,  percentiles)
decomprate56ABBI <- decomprate56ABBI + decomprate45ALL


# 6 -> 7 = dependent on snag fall

#### Why no transition from 6 to 7 in dataset?


# 7 -> 8 = dependent on snag fall
rate78PIEN <- rate56PIEN
decomprate78PIEN <- decompRate(rate78PIEN,  percentiles)
decomprate78PIEN <- decomprate78PIEN + decomprate56PIEN

rate78ABBI <- rate56ABBI
decomprate78ABBI <- decompRate(rate78ABBI,  percentiles)
decomprate78ABBI <- decomprate78ABBI + decomprate56PIEN






# load in dead tree data

deadTreeData <- read.csv("dead_tree_data.csv")
deadTreeData <- deadTreeData %>% 
  select(Transect:Decay)

### Conclass function


conclass <- function(status, decay){
  ifelse(status == 2 & decay %in%  1:2, 3, 
  ifelse(status == 2 & decay %in%  3:4, 4, 
  ifelse(status == 2 & decay %in%  5:6, 5, 
  ifelse(status == 3, 6, 
  ifelse(status == 4 & decay %in%  1:2, 3,
  ifelse(status == 4 & decay %in%  3:4, 4, 
  ifelse(status == 4 & decay %in%  5:6, 5, 
  ifelse(status == 4 & decay == 7, 8,
  ifelse(status == 6, 8, 
  ifelse(status == 5, 7, NA))))))))))
}



deadTreeData$Conclass <- conclass(deadTreeData$Status, deadTreeData$Decay)


decomprates <- rbind(
  decomprate34ALL,
  decomprate45ALL,
  decomprate56PIEN,
  decomprate56ABBI,
  decomprate78PIEN,
  decomprate78ABBI
)


colnames(decomprates) <- percentiles



### merge() here

deadTreeData$DeathYear 







save.image(file = "liveTrees.RData")

