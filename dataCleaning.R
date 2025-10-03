# Data Cleaning

library(dplyr)



# Assumptions made by data:
# (1) 1 dataset with tree identification (plot/tree #)
# (2) DBH, Status, Decay Class, Species, Age columns
# (3) Species column should be cleaned



##### Alpine Dataset
liveTreeData <- read.csv("tree_data.csv")
deadTreeData <- read.csv("dead_tree_data.csv")


liveTreeData <- liveTreeData %>% 
  select(Transect..:Decay, DBH, -Ht.) %>% 
  rename(Transect = Transect..,
         Plot = Plot..,
         Tree = Tree..)
liveTreeData$Status <- as.numeric(liveTreeData$Status)
liveTreeData$DBH <- as.numeric(liveTreeData$DBH)



deadTreeData <- deadTreeData %>% 
  select(Transect:Decay, X2013.DBH, -treeline.pos., -Ht.) %>% 
  rename(Plot = Plot..,
         Tree = Tree..,
         DBH = X2013.DBH)
deadTreeData$Status <- as.numeric(deadTreeData$Status)
liveTreeData$DBH <- as.numeric(liveTreeData$DBH)

fullTreeData <- full_join(liveTreeData, deadTreeData)

fullTreeData$Species <- trimws(fullTreeData$Species)
fullTreeData$Species <- toupper(fullTreeData$Species)

unique(fullTreeData$Species)

fullTreeData <- fullTreeData %>% 
  filter(Species == "ABBI" | Species == "PIEN")

unique(fullTreeData$Species)


write.csv(fullTreeData, file = "AlpineFullTreeData.csv")


##### GUMO Dataset
GUMOlive <- read.csv("GUMOliveTrees.csv")
GUMOdead <- read.csv("GUMOdeadTrees.csv")


GUMOlive <- GUMOlive %>% 
  select(Plot.no.:known.age, -X2004.size.bin) %>% 
  rename(DBH = DBH..2004.,
         Age = known.age)
GUMOlive$Status <- 1
GUMOlive$DecayClass <- NA
GUMOlive$Age[GUMOlive$Age <= 0] <- NA




GUMOdead <- GUMOdead %>% 
  select(Plot.no.:Decay.Class) %>% 
  rename(DBH = Dbh,
         DecayClass = Decay.Class)
GUMOdead$Age <- NA


GUMOfull <- full_join(GUMOlive, GUMOdead)

unique(GUMOfull$Species)

#rename columns to test function
GUMOfull <- GUMOfull %>% 
  filter(Species != "") %>% 
  rename(sp = Species,
         d = DBH,
         a = Age,
         st = Status,
         dc = DecayClass)

unique(GUMOfull$sp)


write.csv(GUMOfull, file = "GumoFullTreeData.csv")



##### Corwina Dataset
corwinaData <- read.csv("corwinaData.csv")
corwinaData <- corwinaData %>% 
  rename(Plot = Plot..) %>% 
  filter(Plot != "C-1")


corwinaAgeData <- read.csv("corwinaAgeData.csv")
corwinaAgeData <- corwinaAgeData %>% 
  select(Plot, Species, DBH, Adjusted.Age)


corwinaData$DBH <- as.numeric(corwinaData$DBH)
corwinaAgeData$DBH <- as.numeric(corwinaAgeData$DBH)



corwinaFullData <- left_join(corwinaData, corwinaAgeData, by = c("Plot", "Species", "DBH"))

corwinaC2Data <- corwinaFullData %>% 
  filter(Plot == "C-2",
         Species %in% c("PIPO", "PSME"))


write.csv(corwinaC2Data, file = "corwinaC2Data.csv")



corwinaCppData <- corwinaFullData %>% 
  filter(Plot != "C-2",
         Species %in% c("PIPO", "PSME"))


write.csv(corwinaCppData, file = "corwinaCppData.csv")



