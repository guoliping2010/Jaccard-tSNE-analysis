library(readxl)
library(openxlsx)
library(stringr)
library(dplyr)
library(proxy)
library(Rtsne)
library(dbscan)
library(ggplot2)

# def Jaccard distance calculation function
calculate_jaccard_distance <- function(data_matrix) {
  jaccard_dist <- as.matrix(dist(data_matrix, method = "Jaccard"))
  diag(jaccard_dist) <- 0 
  return(jaccard_dist)
}

jaccard_tsne_analysis <- function(data_matrix, perplexity = 30, max_iter = 1000, seed=42) {
  jaccard_dist <- calculate_jaccard_distance(data_matrix)
  set.seed(seed)
  tsne_results <- Rtsne(
    X = jaccard_dist,
    is_distance = TRUE,
    perplexity = perplexity,
    max_iter = max_iter,
    verbose = TRUE
  )
  return(tsne_results)
}

#for tsne plot visualization-------
plot_tsne_results <- function(tsne_results, labels = NULL, 
                              title = "t-SNE visualization based on Jaccard distance") {
  plot_data <- data.frame(
    tSNE1 = tsne_results$Y[, 1],
    tSNE2 = tsne_results$Y[, 2]
  )
  
  if (!is.null(labels)) {
    plot_data$Group <- as.factor(labels)
  }
  
  p <- ggplot(plot_data, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = if (exists("Group")) Group else NULL), size = 2) +
    labs(title = title, x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      legend.position = "right"
    )
  
  if (exists("Group", where = plot_data)) {
    p <- p + scale_color_discrete(name = "Group")
  }
  
  print(p)
}


#calculation--------------------------------------------------------
Gene <- read.xlsx("Supplementary_table_8.xlsx", detectDates = TRUE)
rownames(Gene) <- Gene$PID
Gene <- Gene[,c("ASXL1","BCOR","EZH2","RUNX1","SF3B1","SRSF2","STAG2","U2AF1","ZRSR2","TP53")]
Gene[is.na(Gene)] <- 0
cols_to_convert <- colnames(Gene)
Gene[cols_to_convert] <- lapply(Gene[cols_to_convert], as.integer)
Gene <- as.matrix(Gene)

tsne_results <- jaccard_tsne_analysis(Gene, perplexity = 10)
plot_tsne_results(tsne_results) 

tsne_out_data <- as.data.frame(tsne_results$Y)
fit_cluster_dbscan <- dbscan(scale(tsne_out_data[, 1:2]), eps = 0.35)
table(fit_cluster_dbscan$cluster)
tsne_out_data$cl_dbscan <- factor(fit_cluster_dbscan$cluster)

ggplot(tsne_out_data, aes(x = V1, y = V2, colour = cl_dbscan)) +
  geom_point(size = 2, alpha = 2) +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_colour_manual(values = c("#543005", "#D9F0D3", "#B3CDE3", "#1F78B4", "#FF7F00", "#FED9A6", "#C2A5CF", "#FDDAEC", "#276419", "#FBB4AE", "#B2182B"))+
  xlab("TSNE-1") + ylab("TSNE-2")


#elbow method to identify optimal eps-------------------------
embeded <- as.matrix(tsne_results$Y)
dim <- ncol(embeded)  
k <- 2 * dim                 
k_dist <- kNNdist(embeded, k = k)  
sorted_k_dist <- sort(k_dist, decreasing = FALSE)

df <- data.frame(idx = seq_along(sorted_k_dist), dist = sorted_k_dist)
ggplot(df, aes(x = idx, y = dist)) +
  geom_point(size = 0.8) +
  geom_line() +
  geom_hline(yintercept = 0.27, linetype = "dashed", color = "#CAB2D6") +
  geom_hline(yintercept = 0.35, linetype = "dashed", color = "#6A3D9A") +
  labs(x = "Sample index (sorted by distance descending)", 
       y = paste0(k, "th-nearest neighbor distance")) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  annotate("text", x = max(df$idx)*0.95, y = 0.27, label = "eps = 0.27", 
           color = "#CAB2D6", hjust = 1, vjust = -0.5) +
  annotate("text", x = max(df$idx)*0.95, y = 0.35, label = "eps = 0.35", 
           color = "#6A3D9A", hjust = 1, vjust = -0.5)



#different perplexity visualization ------------------------------------------------------------
perplexity_values <- c(5, 10, 20, 30, 50) 
tsne_list <- list()
original_clusters <- tsne_out_data$cl_dbscan
cluster_factor <- factor(original_clusters, levels = sort(unique(original_clusters)))

for (perp in perplexity_values) {
  set.seed(42)
  tsne_out <- jaccard_tsne_analysis(Gene, perplexity = perp)
  coords <- as.data.frame(tsne_out$Y)
  colnames(coords) <- c("tSNE1", "tSNE2")
  coords$cluster <- cluster_factor  
  coords$perplexity <- paste("perp =", perp)   
  
  tsne_list[[as.character(perp)]] <- coords
}


all_tsne <- do.call(rbind, tsne_list)

all_tsne$perplexity <- factor(all_tsne$perplexity, 
                              levels = paste("perp =", perplexity_values))
colors <- c("#543005", "#D9F0D3", "#B3CDE3", "#1F78B4", "#FF7F00", "#FED9A6", "#C2A5CF", "#FDDAEC", "#276419", "#FBB4AE", "#B2182B")
ggplot(all_tsne, aes(x = tSNE1, y = tSNE2, colour = cluster)) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ perplexity, scales = "free", ncol=5) +  
  scale_colour_manual(values = colors, name = "Cluster") +
  theme_minimal() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "TSNE-1", y = "TSNE-2")


