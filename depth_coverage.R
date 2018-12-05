######################################
###Depth of coverage per individual###
######################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after denovo_map.pl, for one analysis

##############################
#Version Stacks 2.08B ou 1.46#
##############################

m <- 3
M <- 2
n <- 1

##from vcf file, option --vcf on populations

setwd("/media/sdb1/RAD/RAD-Requin/stacks.denovo/output_run1_1/")
popinfo <- read.csv("/media/sdb1/RAD/RAD-Requin/info/popmaprun1_65_wp.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])


#Stacks 1.46
recup_vcf <- grep("batch_2.vcf", dir())

res2 <- read.table(dir()[recup_vcf])

for (l in 1:length(sample_names)){
  vcf1 <- res2[, c(9+l)]
  interm <- as.data.frame(t(matrix(unlist(strsplit(as.vector(vcf1), ":")), nrow = 3, ncol = length(vcf1))))
  interm2 <- as.numeric(as.vector(interm[, 2]))
  
  png(sample_names[l], height = 15, width = 20, units = "cm", res = 400)
  plot(table(interm2)[2:length(table(interm2))], main = sample_names[l], type = "l",
       xlab = "coverage", ylab = "nb_loci")
  dev.off()
}

#Stacks 2.08B
recup_vcf <- grep("snps.vcf", dir())

res2 <- read.table(dir()[recup_vcf])

for (l in 1:length(sample_names)){
  vcf1 <- as.vector(res2[, c(9+l)])
  vcf1[which(vcf1 == "./.")] <- ".:.:.:.:."
  interm <- as.data.frame(t(matrix(unlist(strsplit(as.vector(vcf1), ":")), nrow = 5, ncol = length(vcf1))))
  interm2 <- as.numeric(as.vector(interm[, 2]))

  png(sample_names[l], height = 15, width = 20, units = "cm", res = 400)
  plot(table(interm2)[2:length(table(interm2))], main = sample_names[l], type = "l",
       xlab = "coverage", ylab = "nb_loci")
  dev.off()
}


#avant populations

chemin <- "/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini"
setwd(chemin)

recup_indivs <- grep(".matches.tsv", dir())
l <- 10

recup_indiv <- recup_indivs[l]
res1 <- read.table(dir()[recup_indiv])

plot(table(res1[, 7]), main = dir()[l])
