# Fill scales for manuscript

# Cluster scale
cluster_colors <- c("No feature" = "white","Singleton" = "black","Cluster 1" = "green","Cluster 2" = "orange","Cluster 3" = "purple","Cluster 4" = "blue","Cluster 5" = "red","Cluster 6" = "pink","Cluster 7" = "yellow","Cluster 8" = "cyan")
cluster_labels <- names(cluster_colors)
cluster_scale <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels,name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=5, title.position = "top", label.position = "right"))
cluster_scale_v <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels,name="Phylogenetics of resistance", guide = guide_legend(order=3,ncol=2, title.position = "top", label.position = "right"))

# Phylogenetic reconstruction cluster
phylogenetic_reconstruction_cluster_names <- c("Best joint reconstruction",paste0("Reconstruction axis ",1:7))
phylogenetic_reconstruction_cluster_numeric <- c("Best joint",1:7)
phylogenetic_reconstruction_scale <- scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set1"),breaks =  phylogenetic_reconstruction_cluster_names,name = "Reconstruction axis",guide_legend(nrow=1))
phylogenetic_reconstruction_numeric_scale <- scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set1"),breaks =  phylogenetic_reconstruction_cluster_numeric,name = "Reconstruction axis",guide_legend(nrow=1))

# Resistance scale
Resistance_scale <- scale_fill_manual(breaks = c(1,0),values=c("black","white"),labels = c("Resistant","Susceptible"),name="Resistance profile",guide = guide_legend(order=1,title.position = "top", label.position = "right",nrow=2),drop = FALSE )

# Factor scale
factor_scale <- scale_fill_manual(breaks = c(1,0),values=c("black","white"),labels = c("Present","Absent"),name="Genotype",guide = guide_legend(order=3,ncol=1,title.position = "top", label.position = "right"))

# Factor genotypes
factor_genotypes <- function(genotypes){
  factor(genotypes,levels = rev(c("AA552","OmpK36_c25t","OmpK36GD","OmpK36TD","OmpK36_promoter_IS","OmpK36_disruptive","PBP_any","RamR","RamA")))
}

factor_genotypes_boxplot <- function(genotypes){
  factor(genotypes,levels = rev(c("Overall","AA552","OmpK36_c25t","OmpK36GD","OmpK36TD","OmpK36_promoter_IS","OmpK36_disruptive","PBP_any","RamR","RamA")))
}


# Recode functions
recode_genotypes <- function(genotypes){
  recode_factor(genotypes ,
                "Overall" = "Overall",
                "AA552" = "AA552 blaKPC-containing plasmid",
                'OmpK36_c25t'='25 cytosine-to-thymine transition in ompK36',
                'OmpK36GD'='GD insertion in loop 3 of ompK36',
                'OmpK36TD'='TD insertion in loop 3 of ompK36',  
                'OmpK36_promoter_IS'='Insertion sequence at promoter of ompK36',
                "OmpK36_disruptive" = "Putative function-altering variant in ompK36",
                'PBP_any'='Penicillin-binding-protein mutant',
                'RamR'='Putative function-altering variant in ramR efflux pump regulator',
                "RamA"='Putative function-altering variant in ramA efflux pump activator') 
}

collapse_PFAV <- function(genotypes){
  gsub("Putative function-altering variant","PFAV",genotypes)
}

recode_genotypes_for_html <- function(genotypes){
  gsub("ompK36","<i>ompK36</i>",genotypes) %>% gsub("ramR","<i>ramR</i>",.) %>% gsub("ramA","<i>ramA</i>",.) %>% gsub("blaKPC","<i>bla</i><sub>KPC</sub>",.) %>% gsub("variant in","variant <br>in",.)
} 

recode_genotypes_for_complex <- function(genotypes){
  out <- paste0("'", genotypes, "'")
  
  # Line break (plotmath)
  out <- gsub("variant in", "variant~'\n'~in", out)
  
  # Italicize genes (plotmath symbols, not strings)
  out <- gsub("ompK36", "italic(ompK36)", out)
  out <- gsub("ramR",   "italic(ramR)",   out)
  out <- gsub("ramA",   "italic(ramA)",   out)
  out <- gsub("blaKPC", "italic(bla)[KPC]", out)
  
  parse(text = out)
} 

# Recode phenotype clusters
recode_phenotype_clusters <- function(clusters){ 
  recode(clusters, "focal_states" ="Best joint reconstruction",
         "1" ="Reconstruction axis 1",
         "2" ="Reconstruction axis 2",
         "3" ="Reconstruction axis 3",
         "4" ="Reconstruction axis 4",
         "5" ="Reconstruction axis 5",
         "6" ="Reconstruction axis 6",
         "7" ="Reconstruction axis 7")   
}