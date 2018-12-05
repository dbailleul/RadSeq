#################################################
###r80, individual number of loci and coverage###
#################################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after finding good parameters

######################
#Version Stacks 2.08B#
######################

#here on sharks

res2 <- read.table("/media/sdb1/RAD/RAD-Requin/stacks.denovo/output_run1_2/populations.haps.vcf")
popinfo <- read.csv("/media/sdb1/RAD/RAD-Requin/info/popmaprun1_59_wp.csv", sep = ",", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

mat <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = 4))

for(l in 1:length(sample_names)){
  vcf1 <- as.vector(res2[, c(9+l)])
  
  mat[1, l] <- length(vcf1)
  mat[2, l] <- length(vcf1) - length(which(vcf1 == "./."))
  mat[3, l] <- length(which(vcf1 == "./."))
}

res3 <- readLines("/media/sdb1/RAD/RAD-Requin/stacks.denovo/output_run1_2/gstacks.distribs")
posibeg <- grep("BEGIN effective_coverages_per_sample", res3)
posiend <- grep("END effective_coverages_per_sample", res3)
res3 <- res3[(posibeg+2):(posiend-1)]

resf <- as.data.frame(t(matrix(unlist(strsplit(res3, "\t")), nrow = 4, ncol = length(res3))))

if (sum(resf[,1] == sample_names) != length(sample_names)){print("WARNING !!!")}

resf <- as.numeric(as.vector(resf[,2]))
mat[4, ] <- resf

colnames(mat) <- sample_names
rownames(mat) <- c("r80_nb_loci", "indiv_nb_loci", "indiv_no_loci", "indiv_cov")
#write.csv(mat, "/media/sdb1/RAD/RAD-Requin/stacks.denovo/output_run1_2/mat_indiv.csv")

#####################
#Version Stacks 1.46#
#####################

#here on raies
#avec option --vcf_haplotypes dans populations

#setwd("/mnt/sda1/RAD/RAD-Raies/stacks.denovo/output")

res4 <- read.table("/mnt/sda1/RAD/RAD-Raies/stacks.denovo/output/batch_1.haplotypes.vcf")
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_R_mod.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

mat <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = 4))

for(l in 1:length(sample_names)){
  vcf1 <- as.vector(res4[, c(9+l)])
  resf <- as.data.frame(t(matrix(unlist(strsplit(vcf1, ":")), nrow = 2, ncol = length(vcf1))))
  resf <- resf[,1]
  
  mat[1, l] <- length(resf)
  mat[2, l] <- length(resf) - length(which(resf == "./."))
  mat[3, l] <- length(which(resf == "./."))
}

res3 <- readLines("/mnt/sda1/RAD/RAD-Raies/stacks.denovo/output/denovo_map.log")
posibeg <- grep("Depths of Coverage for Processed Samples:", res3)
posiend <- posibeg+length(sample_names)
res3 <- res3[(posibeg+1):(posiend)]

resf <- as.data.frame(t(matrix(unlist(strsplit(res3, ":")), nrow = 2, ncol = length(res3))))

if (sum(resf[,1] == sample_names) != length(sample_names)){print("WARNING !!!")}

resf <- as.data.frame(t(matrix(unlist(strsplit(as.vector(resf[,2]), "x")), nrow = 1, ncol = length(res3))))
mat[4, ] <- as.numeric(as.vector(resf[,1]))

colnames(mat) <- sample_names
rownames(mat) <- c("r80_nb_loci", "indiv_nb_loci", "indiv_no_loci", "indiv_cov")

mat_haps <- mat
#write.csv(mat_haps, "mat_indiv.csv")

