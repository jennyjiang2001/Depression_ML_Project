# Code adapted from online example
# https://osf.io/p9t3y
# https://www.sciencedirect.com/science/article/pii/S1389945722012540

install.packages("summarytools")     
install.packages("networktools")     
install.packages("bootnet")          
install.packages("psychonetrics")   
install.packages("fastshap")         
install.packages("sna")              
install.packages("ltm")              
install.packages("EGAnet")           
install.packages("vip")              
install.packages("huge")             

library(caret)          
library(summarytools)   
library(networktools)   
library(bootnet)        
library(dplyr)          
library(qgraph)         
library(psychonetrics) 
library(fastshap)       
library(sna)            
library(ltm)            
library(EGAnet)         
library(psych)          
library(GPArotation)    
library(vip)          

phqs <- read.csv("final_merged.csv", header = TRUE)  
phqs <- phqs[, 55:63] # read in phq 1-9
print(phqs)                                         

# prepare variable matrix
network_variables <- cbind(phqs)

# read item labels from text file
phq_item <- scan(file = "phq.txt", what = "character", sep = "\n") 

# split dataset into training and testing (70% training)
sample_size <- floor(0.7 * nrow(network_variables))  # compute training set size
set.seed(42)                                         
picked <- sample(seq_len(nrow(network_variables)), size = sample_size)

phq_train <- network_variables[picked, ]             # training set
phq_test  <- network_variables[-picked, ]            # test set


# Estimate network (PMRF) on training set
network <- estimateNetwork(phq_train,
                           default = "ggmModSelect",     
                           corMethod = "npn",            
                           tuning = 0.5,                 
                           labels = phq_item,            
                           missing = "listwise",         
                           verbose = FALSE)              


# Create adjacency matrix from network
adjacency <- 1 * (network$graph != 0)


# Build graphical model using adjacency from training
model_train <- ggm(phq_train, omega = adjacency)
model_train 

# Confirm model structure on testing data
confirmatory <- ggm(phq_test, omega = adjacency)
confirmatory


# Re-estimate network on test data
network2 <- estimateNetwork(phq_test,
                            default = "ggmModSelect",
                            corMethod = "npn",
                            tuning = 0.5,
                            labels = phq_item,
                            missing = "listwise",
                            verbose = FALSE)
adjacency2 <- 1 * (network2$graph != 0)
model_test <- ggm(phq_test, omega = adjacency2)

# Compare correlation structures between train/test
mat1 <- as.data.frame(network$graph)
mat2 <- as.data.frame(network2$graph)

cor_mat1 <- mat1[upper.tri(mat1, diag = FALSE)]  
cor_mat2 <- mat2[upper.tri(mat2, diag = FALSE)]

cor(cor_mat1, cor_mat2, method = "spearman")  

# Visualize network graph with proper labels
# phq_labels <- c("PHQ-1", "PHQ-2", "PHQ-3", "PHQ-4", "PHQ-5", "PHQ-6", "PHQ-7", "PHQ-8", "PHQ-9")

pdf("Network.pdf")  
phq_graph <- qgraph(network$graph,
                    labels = phq_labels,
                    legend = FALSE,
                    vsize = 6.5,
                    label.scale = FALSE,
                    label.cex = 0.8,
                    theme = "colorblind")
dev.off()


# Extract centrality indices (e.g., Expected Influence)
CI <- centrality(phq_graph)


# Plot Expected Influence centrality
pdf("CI.pdf", height = 4, width = 4)
centralityPlot(phq_graph,
               include = "ExpectedInfluence",
               orderBy = "ExpectedInfluence",
               decreasing = FALSE)
dev.off()
