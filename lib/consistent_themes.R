# Fill scales for manuscript

# Cluster scale
cluster_colors <- c("No feature" = "white","Singleton" = "black","Cluster 1" = "green","Cluster 2" = "orange","Cluster 3" = "purple","Cluster 4" = "blue","Cluster 5" = "red","Cluster 6" = "pink","Cluster 7" = "yellow","Cluster 8" = "cyan")
cluster_labels <- names(cluster_colors)
cluster_scale <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels,name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=5, title.position = "top", label.position = "right"))
cluster_scale_v <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels,name="Phylogenetics of resistance", guide = guide_legend(order=3,ncol=2, title.position = "top", label.position = "right"))

# Phylogenetic reconstruction cluster
phylogenetic_reconstruction_cluster_names <- c("Best joint reconstruction",paste0("Reconstruction axis ",1:5))
phylogenetic_reconstruction_cluster_numeric <- c("Best joint",1:5)
phylogenetic_reconstruction_scale <- scale_fill_manual(values = RColorBrewer::brewer.pal(6,"Set1"),breaks =  phylogenetic_reconstruction_cluster_names,name = "Phenotype axis",guide_legend(nrow=1))
phylogenetic_reconstruction_numeric_scale <- scale_fill_manual(values = RColorBrewer::brewer.pal(6,"Set1"),breaks =  phylogenetic_reconstruction_cluster_numeric,name = "Phenotype axis",guide_legend(nrow=1))

# Resistance scale
Resistance_scale <- scale_fill_manual(breaks = c(1,0),values=c("black","white"),labels = c("Resistant","Susceptible"),name="Resistance profile",guide = guide_legend(order=1,title.position = "top", label.position = "right",nrow=2),drop = FALSE )

# Factor scale
factor_scale <- scale_fill_manual(breaks = c(1,0),values=c("black","white"),labels = c("Present","Absent"),name="Genotype",guide = guide_legend(order=3,ncol=1,title.position = "top", label.position = "right"))