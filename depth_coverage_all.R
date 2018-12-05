######################################
###Depth of coverage per individual###
######################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after findings good parameters

##############################
#Version Stacks 2.08B ou 1.46#
##############################

m <- 3
M <- 3
n <- 3

##from vcf file, option --vcf on populations

setwd("/mnt/sda1/RAD/RAD-Raies/stacks.denovo/output")

#Stacks 1.46
recup_vcf <- grep("batch_1.vcf", dir())

#Stacks 2.08B
recup_vcf <- grep("snps.vcf", dir())

res2 <- read.table(dir()[recup_vcf])
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_R_mod.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

for (l in 1:length(sample_names)){
  vcf1 <- res2[, c(9+l)]
  interm <- as.data.frame(t(matrix(unlist(strsplit(as.vector(vcf1), ":")), nrow = 3, ncol = length(vcf1))))
  interm2 <- as.numeric(as.vector(interm[,2]))

  png(sample_names[l], height=15, width=20, units="cm", res=400)
  plot(table(interm2)[2:length(table(interm2))], main = sample_names[l], type = "l",
       xlab = "coverage", ylab = "nb_loci")
  dev.off()
}


#avant populations

recup_indivs <- grep(".matches.tsv", dir())
l <- 10

recup_indiv <- recup_indivs[l]
res1 <- read.table(dir()[recup_indiv])

plot(table(res1[, 7]), main = dir()[l])
