GumoTreeData <- read.csv("GumoFullTreeData.csv")




b <- ageLM(GumoTreeData, speciesCol = sp, ageCol = a, DBHcol = d, measYear = 2004)

b2 <- refLiveDBH(b, refYear = 1922, measYear = 2004, DBHcol = d, speciesCol = sp, avgIncVec = c(PIED = 0.102439024390244/5, PIPO = 0.207317073170732/5, PIST = 0.136585366/5, PSME = 0.190243902439024/5, QUGA = 0.075609756097561/5))

b3 <- conclass(b2, st, dc)

b4 <- decompRate(b3, speciesCol = sp, DBHcol = d, measYear = 2004, refYear = 1922, avgIncVec = c(PIED = 0.102439024390244/5, PIPO = 0.207317073170732/5, PIST = 0.136585366/5, PSME = 0.190243902439024/5, QUGA = 0.075609756097561/5))


