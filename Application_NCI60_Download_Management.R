# data_management:

# go to webpage: https://discover.nci.nih.gov/cellminer/loadDownload.do

# download to "folderapplication" the two data sets:

# "RNA: Affy HG-U133(A-B)" with GCRMA
# "Protein: Antibody Array DTB" with log2




#setwd("//folderapplication")


library(readxl)
NCI60 <- read_excel("RNA__Affy_HG_U133(A_B)_GCRMA.xls")
KRT19_protein <- read_excel("Protein__Antibody_Array_DTB_log2.xls")
number_cancer_cell_lines <- 59


# predictors:
NCI_60_pred <- NCI60[-c(1:10),]
colnames(NCI_60_pred) <- NCI60[10,]

# delete column with empty entries:

NCI_60_pred <- NCI_60_pred[,-which(colnames(NCI_60_pred)=="LC:NCI-H23")]

# log2:
for (i in 8:66) {
  NCI_60_pred[[i]] <- as.numeric(NCI_60_pred[[i]])
}
NCI_60_pred <- as.data.frame(NCI_60_pred)

NCI_60_pred$`Gene name d`[which(NCI_60_pred$`Gene name d`=="-")] <- paste("ident_",NCI_60_pred$`Identifier c`[which(NCI_60_pred$`Gene name d`=="-")],sep="")

#which(NCI_60_income$`Gene name d`=="E2F4")

# delete column with negative entries:
NCI_60_pred <- NCI_60_pred[-16206,]
for (i in 8:66) {
  NCI_60_pred[,i] <-NCI_60_pred[,i]
}



uniname <- unique(NCI_60_pred$`Gene name d`)

pred_matrix <- matrix(0,ncol=number_cancer_cell_lines,nrow = length(uniname))

# median of multiple studies for gene expression data with same name
for (j in 8:66) {
  for (i in 1:length(uniname)) { pred_matrix[i,j-7] <- median(c(NCI_60_pred[which(NCI_60_pred$`Gene name d`==uniname[i]),j]))
  }
}

trans_pred_matrix <- t(pred_matrix)



# outcome:
#which(KRT19_protein[,1]=="KRT19")
KRT19_protein_vector <- as.vector(KRT19_protein[which(KRT19_protein[,1]=="KRT19"),-c(1:6)])

KRT19_protein_vector <- as.numeric(KRT19_protein_vector)

# delete "LC:NCI-H23" entry:

co_names<- as.character(KRT19_protein[10,7:66])
co_names2 <- co_names[-c(40)]

KRT19_protein_vector <- KRT19_protein_vector[-c(40)]

KRT19_protein_vector2 <- 2^KRT19_protein_vector# getting the original size without log2 transformation

#matching:
#co_names(NCI_60_pred)
#co_names2 # already same position

# names are similar and same order
data_set <- data.frame(KRT19_protein_vector2,trans_pred_matrix)



colnames(data_set) <- c("KRT19_protein",uniname)
colnames(data_set) <- gsub("-", "_", colnames(data_set))


#setwd("//folderapplication")
save(data_set,file="data_set_application.RData")