#################################################
###r80, individual number of loci and coverage###
#################################################

#Diane, on March 2018
#diane.bailleul.pro@gmail.com

#after setting the goos parameters

##############################
#Version Stacks 2.08B ou 1.46#
##############################

chemin <- c("/mnt/sda1/RAD/RAD-Raies/stacks.denovo/output")
setwd(chemin)

popinfo <- read.csv("/mnt/sda1/RAD/RAD-Raies/info/popmaprun2_R_mod.csv", sep = "\t", header = FALSE)
sample_names <- as.vector(popinfo[, 1])

n <- 3
nindiv <- length(sample_names)
mat_recup <- as.data.frame(matrix(0, ncol = nindiv, nrow = 1))
m <- 3
M <- 3

snps <- read.table(gzfile("batch_1.catalog.snps.tsv.gz"))
tags <- read.table(gzfile("batch_1.catalog.tags.tsv.gz"))

mat <- as.data.frame(matrix(0, ncol = nindiv, nrow = length(unique(snps[, 3]))))
rownames(mat) <- paste("A_", rownames(mat), sep = "")

k <- 1
for (i in 1:nrow(snps)){
  if (i == 1){
    recup <- which(tags[, 3] == snps[i, 3])
    rownames(mat)[k] <- snps[i, 3]
    recup2 <- unlist(strsplit(as.vector(tags[recup, 8]), ","))
    for (j in 1:nindiv){
      indice <- paste(j, "_", sep = "")
      mat[k, j] <- length(grep(indice, recup2))
    }
    k <- k+1
    
  } else {
    
    recup <- which(tags[, 3] == snps[i, 3])
    if (snps[i, 3] != snps[i-1, 3]){
      rownames(mat)[k] <- snps[i, 3]
      recup2 <- unlist(strsplit(as.vector(tags[recup, 8]), ","))
      for (j in 1:nindiv){
        indice <- paste(j, "_", sep = "")
        mat[k, j] <- length(grep(indice, recup2))
      }
      k <- k+1
    }
    
  }
}

mat_loci <- mat

for (l in 1:nindiv){
  
  mat_indiv <- as.data.frame(matrix(0, ncol = nindiv, nrow = nrow(mat_loci)))
  rownames(mat_indiv) <- rownames(mat_loci)
  
  for (o in 1:nrow(mat_loci)){
    loc <- mat_loci[o,]
    if (loc[l] != 0){
      recup <- sum(loc != 0)
      #mat_indiv[o, recup] <- 1
      mat_indiv[o, recup] <- loc[l]
    }
  }
  plot(apply(mat_indiv, 2, sum), type = "l", col = "red", main = paste("m=", m, " M=", M, " n=", n, " indiv=", l, sep = ""))
  mat_recup <- rbind(mat_recup, apply(mat_indiv, 2, sum))
  rownames(mat_recup)[nrow(mat_recup)] <- paste("m", m, "_M", M, "_n", n, "_l", l, sep = "") 
}

mat_recup <- mat_recup[-1,]

#write.csv(mat_recup, "overlapping_loci_details.csv")


for (l in 1:nindiv){

  png(paste("m=", m, " M=", M, " n=", n, " indiv=", sample_names[l], sep = ""), height = 15, width = 20, units = "cm", res = 400)
  plot(as.numeric(mat_recup[l,]), type = "l", col = "red", main = paste("m=", m, " M=", M, " n=", n, " indiv=", sample_names[l], sep = ""))
  dev.off() 

}