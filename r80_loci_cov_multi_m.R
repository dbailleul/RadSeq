#################################################
###r80, individual number of loci and coverage###
#################################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after denovo_map.pl, for several analysis to choose m

#####################
#Version Stacks 1.46#
#####################

#here on raies
#avec option --vcf_haplotypes dans populations

chemin <- "/mnt/sda1/RAD/RAD-Raies/tests.denovo/"
setwd(chemin)
popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_Rmini_mod.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

m_min <- 2
m_max <- 8
n <- 1
M <- 2

folder_name <- "output_mini_m"

mat1 <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = length(m_min:m_max)))
mat2 <- as.data.frame(matrix(0, ncol = length(m_min:m_max), nrow = 1))
mat3 <- as.data.frame(matrix(0, ncol = length(sample_names), nrow = length(m_min:m_max)))


for (m in m_min:m_max){

  res4 <- read.table(paste(chemin, folder_name, m, "/batch_2.haplotypes.vcf", sep = ""))

  for (l in 1:length(sample_names)){
    vcf1 <- as.vector(res4[, c(9+l)])
    resf <- as.data.frame(t(matrix(unlist(strsplit(vcf1, ":")), nrow = 2, ncol = length(vcf1))))
    resf <- resf[,1]
  
    mat1[which(m_min:m_max == m), l] <- length(resf) - length(which(resf == "./."))

  }

mat2[1, which(m_min:m_max == m)] <- length(resf)

}

rownames(mat1) <- m_min:m_max
colnames(mat1) <- sample_names

colnames(mat2) <- m_min:m_max
rownames(mat2) <- "r80"

#write.csv(mat1, "nb_loci_m_variations.csv")
#write.csv(mat2, "r80_haps_m_variations.csv")


mat4 <- as.data.frame(matrix(0, ncol = m_min:m_max, nrow = 2))

for (m in m_min:m_max){
  
  res3 <- readLines(paste(chemin, folder_name, m, "/batch_2.populations.log", sep = ""))
  posibeg <- grep("# Distribution of the number of SNPs per locus.", res3)
  posiend <- length(res3)
  res3 <- res3[(posibeg+2):(posiend)]
  
  resf <- as.data.frame(t(matrix(unlist(strsplit(res3, "\t")), nrow = 2, ncol = length(res3))))
  resf <- as.numeric(as.vector(resf[,2]))
  
  mat4[1, which(m_min:m_max == m)] <- sum(resf)
  mat4[2, which(m_min:m_max == m)] <- sum(resf[-1])
  
}

colnames(mat4) <- m_min:m_max
rownames(mat4) <- c("all_r80", "poly_r80")

#write.csv(mat4, "r80_batch_m_variations.csv")


for (m in m_min:m_max){

  res3 <- readLines(paste(chemin, folder_name, m, "/denovo_map.log", sep = ""))
  posibeg <- grep("Depths of Coverage for Processed Samples:", res3)
  posiend <- posibeg+length(sample_names)
  res3 <- res3[(posibeg+1):(posiend)]

  resf <- as.data.frame(t(matrix(unlist(strsplit(res3, ":")), nrow = 2, ncol = length(res3))))

    if (sum(resf[,1] == sample_names) != length(sample_names)){print("WARNING !!!")}

  resf <- as.data.frame(t(matrix(unlist(strsplit(as.vector(resf[,2]), "x")), nrow = 1, ncol = length(res3))))

  mat3[which(m_min:m_max == m), ] <- as.numeric(as.vector(resf[,1]))

}

rownames(mat3) <- m_min:m_max
colnames(mat3) <- sample_names

#write.csv(mat3, "cov_m_variations.csv")