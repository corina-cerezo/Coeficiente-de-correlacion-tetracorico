library(psych)
?tetrachoric
tetrachoric(matrix(c(44268,193,14,0),2,2))  #MPLUS reports.24

tetrachoric(matrix(c(1562,42,383,94),2,2))  #Ejemplo de Ekrtöm
rt <- tetrachoric(matrix(c(631,125,147,147),2,2))  #ilustración 1 de Pearson
rt
rt$rho
rt$tau
chisq.test(matrix(c(1562,42,383,94),2,2))



rt <- tetrachoric(matrix(c(1766,842,842,722),2,2))  #ilustración 3 de Pearson
rt$rho
rt$tau
rt <- tetrachoric(matrix(c(254,136,156,193),2,2))  #ilustración 4 de Pearson
rt
rt$rho
rt$tau
rt$objective
rt$Call
data(bock)
responses <- table2df(bock.table[,2:6],count=bock.table[,7],
                      labs= paste("lsat6.",1:5,sep=""))
describe(responses)
