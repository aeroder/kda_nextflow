#You can either submit a large as a KDA input or can break this list up by Membership to run in parallel if speed is an issue. 
KDA_DF = data.frame(Gene = NA, Membership = NA)

##################
#WGCNA Modules
##################
Dataset = "Colorectal-Selected-Meds_23Q1"
# Dataset = "mCRPC_22Q4"
# Dataset = "NSCLC-selected-meds_22Q4"
# Dataset = "Ovarian_23Q1"
# Dataset = "Pancreas-Selected-Meds_23Q1"

#WGCNA Output to KDA input. 
OuputFolder = paste0("/Users/phillipcomella/Library/CloudStorage/OneDrive-Pathos/p0044_coexpression_nets_workflow/results/", Dataset, "/tables/")
OuputFolderKDA = paste0(OuputFolder, "KDA")
if (!dir.exists(OuputFolderKDA)){
  dir.create(OuputFolderKDA)
}

#This needs to match whatever the Bayesian Network is
hgnc_symbolORensembl_gene_id = "ensembl_gene_id"
WGCNAgeneInfoFile = read.csv(paste0(OuputFolder, "wgcna_gene_info.csv"))
WGCNAgeneInfoFile = WGCNAgeneInfoFile[,c(hgnc_symbolORensembl_gene_id, "module_color")]
colnames(WGCNAgeneInfoFile) = c("Gene", "Membership")
KDA_DF = rbind(KDA_DF, WGCNAgeneInfoFile)

##################
#Signatures
##################

#os_signatures
os_signatures = readRDS(paste0(OuputFolder,"os_signatures.Rds"))
hist((os_signatures$statistic))
qunatile = quantile(abs(os_signatures$statistic), probs = seq(0, 1, 0.05))

#Only Significant genes
os_signatures = os_signatures[which(os_signatures$adj.p <=0.05 & abs(os_signatures$statistic)>=qunatile[20]),]

#Separate based on Direction
os_signatures_UP = os_signatures[os_signatures$statistic>0,]
os_signatures_UP$Membership = "os_signatures_UP"
os_signatures_UP = os_signatures_UP[,c("ensembl_gene_id","Membership")]
colnames(os_signatures_UP) = c("Gene", "Membership")
KDA_DF = rbind(KDA_DF, os_signatures_UP)

os_signatures_DOWN = os_signatures[os_signatures$statistic<0,]
os_signatures_DOWN$Membership = "os_signatures_DOWN"
os_signatures_DOWN = os_signatures_DOWN[,c("ensembl_gene_id","Membership")]
colnames(os_signatures_DOWN) = c("Gene", "Membership")
KDA_DF = rbind(KDA_DF, os_signatures_DOWN)


#diff_expression
#Many of these tables seem very large, pausing to use them until we can find a good way to reduce them. 
# diff_expression = readRDS(paste0(OuputFolder,"diff_expression.Rds"))
# list =5
# for(list in 1:length(diff_expression)){
#   x = as.data.frame(diff_expression[[list]])
#   hist(x$adj.p.value)
#   if(nrow(x)>1){
#     x$ensembl_gene_id
#   }
# }

##################
#Global. Technically this should probably be done in the KDA R code as it just needs the BN to call the entire BN but this will also work because the WGCNA input file contains the entire dataset and the BN typically is just a subset of the entire dataset. 
##################
global = data.frame(Gene = unique(KDA_DF$Gene), Membership = "Global")
KDA_DF = rbind(KDA_DF, global)


KDA_DF = na.omit(KDA_DF)
write.table(KDA_DF,paste0(OuputFolderKDA, "/KDAInputFile.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
