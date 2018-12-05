#################################
###Filter Quality on Raw Reads###
#################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after demultiplex, after merging lines

#here on BSH run 1 & 2

data_run1 <- read.csv("run1_qual.csv")
head(data_run1)

data_run1 <- data_run1[order(data_run1[,3]),]

plot(xlim = c(100000, 17000000), x = data_run1[,3], y = 1:95, xlab = "Read_Number",
     yaxt = "n", ylab = "Individuals",
     pch = 20, col = data_run1[,2],
     xaxt = "n")
axis(side = 2, at = 1:95, labels = data_run1[,1], las = 1, cex.axis = 0.5)
axis(side = 1, at = c(100000,1000000,5000000,10000000,15000000, 20000000), las=2, cex.axis = 0.8)

legend(x = 12000000, y = 40,
       title = "Populations",
       legend = unique(data_run1[,2]), fill = unique(data_run1[,2]))

data_run2 <- read.csv("run2_qual.csv")
head(data_run2)

data_run2 <- data_run2[order(data_run2[,3]),]
#de 600,000 à 30,000,000

plot(xlim = c(600000, 30000000), x = data_run2[,3], y = 1:96, xlab = "Read_Number",
     yaxt = "n", ylab = "Individuals",
     pch = 20, col = data_run2[,2],
     xaxt = "n")
axis(side = 2, at = 1:96, labels = data_run2[,1], las = 1, cex.axis = 0.5)
axis(side = 1, at = c(600000,1000000,5000000,10000000,15000000, 20000000, 30000000), las=2, cex.axis = 0.8)

legend(x = 20000000, y = 40,
       title = "Populations",
       legend = unique(data_run2[,2]), fill = unique(data_run2[,2]))