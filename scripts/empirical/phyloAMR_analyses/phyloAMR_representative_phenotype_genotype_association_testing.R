setwd('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Phyloaware/joint-uncertainty-ms/')

#### ENVIRONMENT ####
#Packages
packages <- c("tidyverse",'corHMM','phyloAMR','pbmcapply') 

#Load Packages
invisible(lapply(packages,library,character.only=T,quietly=T))

# Source functions
source("./scripts/empirical/phyloAMR_analyses/phyloAMR_joint_uncertainty_transition_association_testing_algorithm.R")

# Print environment
set.seed(0104) 
Sys.info()
sessionInfo()

#### DATA #####

dat <- readRDS("./data/empirical/dataset/df.RDS")
phy <- read.tree("./data/empirical/tree/tree.treefile") 
genotypes <- readRDS("./data/empirical/resistance_genotypes/resistance_genotypes.RDS") 
representative_parent_child_dfs <- readRDS('./data/empirical/phyloAMR/representative_blbli_joint_uncertainty_parent_child_dfs.RDS') 

# Curate data
dat <- left_join(dat, genotypes)
comparitors <- colnames(genotypes) %>% subset(.!= 'isolate_no')

#### ASSOCIATION TESTING #### 
transition_association_testing <- joint_uncertainty_transition_association_testing(comparitors = comparitors, df = dat, tr = phy, tip_name_variable = 'isolate_no',trait = 'blbli_dich_num', joint_unc_trait_parent_child_dfs = representative_parent_child_dfs, test = 'both', num_permutations = 1000, num_cores = 13)

saveRDS(transition_association_testing,'./data/empirical/phyloAMR/phenotype_genotype_transition_association_statistics.RDS')
