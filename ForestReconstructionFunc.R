# Function for Forest Stand Reconstruction
# Assumes Cleaned Dataset (fullTreeData)


AlpineTreeData <- read.csv("AlpineFullTreeData.csv")
GumoTreeData <- read.csv("GumoFullTreeData.csv")



# Function to run linear regression on aged trees

ageLM <- function(data, Species = Species, Age = Age, DBH = DBH, Year) {
  
  # column names as strings
  Species <- deparse(substitute(Species))
  Age <- deparse(substitute(Age))
  DBH <- deparse(substitute(DBH))
  
  # new dataframe for predicted ages
  ageData <- data
  
  for (sp in unique(data[[Species]])) {
    formula <- as.formula(paste(Age, "~", DBH))
    
    # dataframe for lm (has the target species with Age not NA)
    speciesData <- ageData[ageData[[Species]] == sp & !is.na(ageData[[Age]]), ]
    
    if (nrow(speciesData) > 1) {
      fit <- lm(formula, data = speciesData)
      
      missingRows <- is.na(ageData[[Age]]) & ageData[[Species]] == sp
      
      if (any(missingRows)) {
        predictedAge <- predict(fit, newdata = ageData[missingRows, ])
        ageData[[Age]][missingRows] <- predictedAge
      }
    } 
  }
  
  ageData[[Age]] <- ceiling(as.numeric(ageData[[Age]]))
  ageData$estabYear <- Year - ageData[[Age]]
  
  
  return(ageData)
}



# a <- ageLM(data = AlpineTreeData, Year = 2013)
# b <- ageLM(data = GumoTreeData, species, age, dbh, Year = 2004)

