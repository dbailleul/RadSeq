#################################################
###r80, individual number of loci and coverage###
#################################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after denovo_map.pl, for one analysis

######################
#Version Stacks 2.08B#
######################

#here on raies

res2 <- read.table("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini/populations.haps.vcf")
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/tests.denovo/info/popmaprun2_Rmini.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

mat <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = 4))

for(l in 1:length(sample_names)){
  vcf1 <- as.vector(res2[, c(9+l)])
  
  mat[1, l] <- length(vcf1)
  mat[2, l] <- length(vcf1) - length(which(vcf1 == "./."))
  mat[3, l] <- length(which(vcf1 == "./."))
}

res3 <- readLines("/media/sdb1/RAD/RAD-Requin/tests.denovo/output_NZ_3/gstacks.distribs")
posibeg <- grep("BEGIN effective_coverages_per_sample", res3)
posiend <- grep("END effective_coverages_per_sample", res3)
res3 <- res3[(posibeg+2):(posiend-1)]

resf <- as.data.frame(t(matrix(unlist(strsplit(res3, "\t")), nrow = 4, ncol = length(res3))))

if (sum(resf[,1] == sample_names) != length(sample_names)){print("WARNING !!!")}

resf <- as.numeric(as.vector(resf[,2]))
mat[4, ] <- resf

colnames(mat) <- sample_names
rownames(mat) <- c("r80_nb_loci", "indiv_nb_loci", "indiv_no_loci", "indiv_cov")
#write.csv(mat, "mat_indiv.csv")


#####################
#Version Stacks 1.46#
#####################

#here on raies
#avec option --vcf_haplotypes dans populations

res4 <- read.table("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini/batch_2.haplotypes.vcf")
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_Rmini.csv", sep = "\t", header = FALSE)
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

res3 <- readLines("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini/denovo_map.log")
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


#old

res2 <- readLines("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini_2/batch_2.haplotypes.tsv")
res2[1]
res2[1] <- "Catalog_ID\tCnt\tPtv_GDL_02.all\tPtv_GDL_04.all\tPtv_OI_15.all\tPtv_OI_16.all\tPtv_med_21.all\tPtv_med_23.all\tPv_23.all\tPv_24.all\tRm_01.all\tRm_01b.all"
writeLines(res2, "/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini_2/batch_2.haplotypes.tsv")

res2 <- read.table("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini_2/batch_2.haplotypes.tsv", header = TRUE)
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_Rmini.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

mat <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = 4))

for(l in 1:length(sample_names)){
  vcf1 <- as.vector(res2[, c(2+l)])
  
  mat[1, l] <- length(vcf1)
  mat[2, l] <- length(vcf1) - length(which(vcf1 == "-"))
  mat[3, l] <- length(which(vcf1 == "-"))
}

res3 <- readLines("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini_2/denovo_map.log")
posibeg <- grep("Depths of Coverage for Processed Samples:", res3)
posiend <- posibeg+length(sample_names)
res3 <- res3[(posibeg+1):(posiend)]

resf <- as.data.frame(t(matrix(unlist(strsplit(res3, ":")), nrow = 2, ncol = length(res3))))

if (sum(resf[,1] == sample_names) != length(sample_names)){print("WARNING !!!")}

resf <- as.data.frame(t(matrix(unlist(strsplit(as.vector(resf[,2]), "x")), nrow = 1, ncol = length(res3))))
mat[4, ] <- as.numeric(as.vector(resf[,1]))

colnames(mat) <- sample_names
rownames(mat) <- c("r80_nb_loci", "indiv_nb_loci", "indiv_no_loci", "indiv_cov")
#write.csv(mat, "mat_indiv.csv")

res5 <- read.table("/mnt/sda1/RAD/RAD-Raies/tests.denovo/output_mini_2/batch_2.vcf")
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_Rmini.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

mat <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = 4))

for(l in 1:length(sample_names)){
  vcf1 <- as.vector(res5[, c(9+l)])
  resf <- as.data.frame(t(matrix(unlist(strsplit(vcf1, ":")), nrow = 3, ncol = length(vcf1))))
  resf <- resf[,1]
  
  mat[1, l] <- length(resf)
  mat[2, l] <- length(resf) - length(which(resf == "./."))
  mat[3, l] <- length(which(resf == "./."))
}
