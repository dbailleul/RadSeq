########################
###r80 for curves M/n###
########################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after denovo_map.pl, for several analysis to choose M/n, given m

#####################
#Version Stacks 1.46#
#####################

#here on raies
#avec option --vcf_haplotypes dans populations

chemin <- "/media/sdb1/RAD/RAD-Requin/tests.denovo/"
setwd(chemin)
popinfo <- read.table("/media/sdb1/RAD/RAD-Requin/info/popmaprun12_21_wp", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

m <- 3

Mn_min <- 1
Mn_max <- 9

folder_name <- "output_opti_m"

mat1 <- as.data.frame(matrix(0, ncol = Mn_min:Mn_max, nrow = 2))

for (Mn in Mn_min:Mn_max){

  res3 <- readLines(paste(chemin, folder_name, m, "_Mn", Mn, "/populations.distribs", sep = ""))
  posibeg <- grep("BEGIN snps_per_loc_postfilters", res3)
  posiend <- grep("END snps_per_loc_postfilters", res3)
  res3 <- res3[(posibeg + 3):(posiend - 1)]

  resf <- as.data.frame(t(matrix(unlist(strsplit(res3, "\t")), nrow = 2, ncol = length(res3))))
  resf <- as.numeric(as.vector(resf[,2]))
  
  mat1[1, which(Mn_min:Mn_max == Mn)] <- sum(resf)
  mat1[2, which(Mn_min:Mn_max == Mn)] <- sum(resf[-1])

}

colnames(mat1) <- Mn_min:Mn_max
rownames(mat1) <- c("all_r80", "poly_r80")

#write.csv(mat1, "r80_batch_m3_Mn_variations.csv")

png(paste("r80_batch_m = ", m, sep = ""), height = 15, width = 20, units = "cm", res = 400)
plot(x = Mn_min:Mn_max, y = mat1[1,], type = "l", ylim = c(0, max(mat1[1,])), col = "royalblue",
    ylab = "r80_loci", main = paste("m = ", m, sep = ""))
lines(x = Mn_min:Mn_max, y = mat1[2,], type = "l", ylim = c(0, max(mat1[2,])), col = "red")
dev.off()

#From haplotypes vcf

mat2 <- as.data.frame(matrix(0, ncol = Mn_min:Mn_max, nrow = 1))

for (Mn in Mn_min:Mn_max){
  
  res4 <- read.table(paste(chemin, folder_name, m, "_Mn", Mn, "/populations.haps.vcf", sep = ""))
  mat2[1, which(Mn_min:Mn_max == Mn)] <- nrow(res4)
  
}

colnames(mat2) <- Mn_min:Mn_max
rownames(mat2) <- "r80"

#write.csv(mat2, "r80_haps_m3_Mn_variations.csv")

png(paste("r80_haps_m = ", m, sep = ""), height = 15, width = 20, units = "cm", res = 400)
plot(x = Mn_min:Mn_max, y = mat2[1,], type = "l", ylim = c(0, max(mat1[1,])), col = "royalblue",
     ylab = "r80_loci", main = paste("m = ", m, sep = ""))
dev.off()
