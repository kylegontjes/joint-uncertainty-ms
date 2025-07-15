joint_uncertainty_transition_association_testing <- function(comparitors, df, tr, tip_name_variable, trait, joint_unc_trait_parent_child_dfs, test='both', num_permutations = 1000, num_cores = 7){
  df <- df[match(tr$tip.label, df[[tip_name_variable]]), ]
  tr$node.label <- NULL
  
  trait_names <- names(joint_unc_trait_parent_child_dfs)
  
  # Trait  data curation
  print(paste0("Starting trait data curation at: ", Sys.time())) 
  
  if(test %in% c("downstream","both")){
    stretches <- lapply(joint_unc_trait_parent_child_dfs,FUN=function(x){
      phyloAMR:::get_trait_traces_on_tree(parent_child_df = x, tr = tr, node_states = 'joint')
    })
    names(stretches) <- trait_names 
  }
  
  print(paste0("Finished trait data curation at: ", Sys.time()))
 
  # Genotype data curation
  print(paste0("Started running corHMM on genotypes: ", Sys.time()))
  genotype_asr_objs <- pbmcapply::pbmclapply(comparitors,FUN=function(x){
    asr_simplified(df = df, tr = tr, tip_name_variable = tip_name_variable, trait = x)
  },mc.cores=num_cores)
  names(genotype_asr_objs) <- comparitors 
  
  gc()
  print(paste0("Finished corHMM on genotypes at: ", Sys.time()))
  
  # Observation data
  ## Synchronous testing
  if(test %in% c("synchronous","both")){
    print(paste0("Started running synchronous on genotypes: ", Sys.time()))

    synchronous_obs_results <- lapply(joint_unc_trait_parent_child_dfs, curate_synchronous_data, genotype_asr_objs=genotype_asr_objs, comparitors = comparitors, num_cores = num_cores) 
    names(synchronous_obs_results) <- trait_names
    print(paste0("Done w/ running synchronous on genotypes: ", Sys.time()))
  }
  
  ## Downstream testing
  if(test %in% c("downstream","both")){
    print(paste0("Started running downstream on genotypes: ", Sys.time()))
    
    downstream_obs_results <- lapply(trait_names, curate_downstream_data, stretches = stretches,genotype_asr_objs = genotype_asr_objs, comparitors = comparitors, num_cores = num_cores)
    names(downstream_obs_results) <- trait_names
    
    print(paste0("Done w/ running downstream on genotypes: ", Sys.time()))
  }
   
  # Generate permutations
  ## Permuted traits
  print(paste0("Started permuting features of genotypes: ", Sys.time()))
  
  num_isolates <- nrow(df)
  permutations <-  pbmclapply(comparitors,FUN=function(x){
    generate_permutations(comparitor = x, df = df , num_permutations = num_permutations, tip_name_variable = tip_name_variable, num_isolates=num_isolates)
  }, mc.cores = num_cores)
  names(permutations) <- comparitors
  
  print(paste0("Done w/ permuting features of genotypes: ", Sys.time()))
  gc()
  
  ## Permuted rates
  print(paste0("Started permuting ancestral states: ", Sys.time()))
  tip_names <- df[[tip_name_variable]]
  permuted_rates <- pbmclapply(comparitors,FUN=function(x){
    generate_permuted_rates(genotype =  x, genotype_asr_objs=genotype_asr_objs, permutations=permutations, df=df, tr=tr, tip_names=tip_names)
  }, mc.cores=num_cores)
  names(permuted_rates) <- comparitors
  
  print(paste0("Done w/ permuting ancestral states: ", Sys.time()))
  
  # Clean before empirical stats
  rm(genotype_asr_objs)
  gc() 
  
  # Generate empirical test statistics
  ## Synchronous
  if(test %in% c("synchronous","both")){
    print(paste0("Started running syncronous pvalue testing: ", Sys.time()))
    
    synchronous_pvalue_test_results <- lapply(trait_names,FUN=function(x){
      pbmclapply(comparitors,FUN=function(y){
        generate_synchronous_test_pvalue(trait_synchronous_obs = synchronous_obs_results[[x]][synchronous_obs_results[[x]][["comparitor"]] == y,], trait_parent_child = joint_unc_trait_parent_child_dfs[[x]],permutation_parent_child =  permuted_rates[[y]], num_permutations = num_permutations)
      },mc.cores = num_cores) %>% `names<-`(comparitors)})
    names(synchronous_pvalue_test_results) <- trait_names 
     
    print(paste0("Done w/ syncronous pvalue testing: ", Sys.time()))
    
  }
  
  ## Downstream
  if(test %in% c("downstream","both")){
    print(paste0("Started running downstream pvalue testing: ", Sys.time()))
    
    downstream_pvalue_test_results <- lapply(trait_names,FUN=function(x){
      pbmclapply(comparitors,FUN=function(y){
        generate_downstream_test_pvalue(trait_downstream_result = downstream_obs_results[[x]][downstream_obs_results[[x]][["comparitor"]] == y,],trait_parent_child = trait_dfs[[x]],permutation_parent_child = permuted_rates[[y]], stretches = stretches[[x]], num_permutations = num_permutations) 
      },mc.cores = num_cores) %>% `names<-`(comparitors)})
    names(downstream_pvalue_test_results) <- trait_names 
    
    print(paste0("Done w/ downstream pvalue testing: ", Sys.time()))
    
  }
   
  rm(permuted_rates)
  gc()
  
  GWAS_results <- list(
    trait=trait,
    comparitors = comparitors,
    node_states='joint',
    num_permutations=num_permutations, 
    synchronous_results = synchronous_obs_results,
    synchronous_permutations = synchronous_pvalue_test_results,
    downstream_results=downstream_obs_results,
    downstream_permutations = downstream_pvalue_test_results
  )
  
  print(paste0("GWAS done at: ", Sys.time()))
  
  return(GWAS_results)
  
}

asr_simplified <- function(df, tr, tip_name_variable, trait, upper_bound = 1e50, lower_bound = 1e-9){
  # Run corHMM 
  print("Performing model finding to use best fitting model")
  
  find_best_model <-  find_best_asr_model(df = df, tr = tr, tip_name_variable =  tip_name_variable, trait =  trait, node_states = 'joint',upper_bound = upper_bound,lower_bound = lower_bound)  
  corHMM_out <- find_best_model[['corHMM_output_best_model']]
  rm(find_best_model)
  
  index.mat <- corHMM_out$index.mat
  solution <- corHMM_out$solution
  states <- corHMM_out$states
  rm(corHMM_out)
  
  #Get parent child
  outcome_str <- df[, trait] %>% `names<-`(df[[tip_name_variable]])
  parent_child_df <- get_parent_child_data(tr = tr, ancestral_states = states, trait_data =  outcome_str, confidence_threshold =  NULL, node_states = 'joint')
  parent_child_df <- get_continuation_data(parent_child_df, node_states = 'joint')
  
  rm(states)
  
  # Return results
  results <- list(parent_child_df = parent_child_df, index.mat = index.mat, solution = solution)
  return(results)
}

curate_synchronous_data <- function(trait_asr,genotype_asr_objs,comparitors,num_cores){
  # Run synchronous testing
  synchronous_obs_data <- pbmclapply(genotype_asr_objs,FUN=function(x){
    synchronous_transitions(comparitor_parent_child_df = x[['parent_child_df']], trait_parent_child_df = trait_asr, node_states = 'joint', confidence = NULL)
  },mc.cores=num_cores)
  # Curate data
  names(synchronous_obs_data) <- comparitors 
  synchronous_obs_data_df <- do.call(rbind.data.frame,synchronous_obs_data)
  synchronous_obs_data_df$comparitor = comparitors 
  # Remove 
  rm(synchronous_obs_data)
  return(synchronous_obs_data_df)
}

curate_downstream_data <- function(cluster,stretches,genotype_asr_objs,comparitors,num_cores){
  # Run downstream testing
  downstream_obs_data <- pbmclapply(genotype_asr_objs,FUN=function(x){
    phyloAMR:::get_gain_loss_on_stretches(comparitor_parent_child_df = x[['parent_child_df']], downstream_nodes = stretches[[cluster]], node_states = 'joint', confidence = NULL)
  },mc.cores=num_cores)  
  # Curate data
  names(downstream_obs_data) <- comparitors 
  downstream_obs_data_df <- do.call(rbind.data.frame,downstream_obs_data)
  downstream_obs_data_df$comparitor = comparitors  
  # Remove
  rm(downstream_obs_data)
   
  return(downstream_obs_data_df)
}


generate_permutations <- function(comparitor,df,num_permutations,tip_name_variable,num_isolates=num_isolates) {
  comparitor_vals <- df[[comparitor]]
  trait_runs <- replicate(num_permutations,comparitor_vals[sample.int(num_isolates)],simplify = F)
  return(trait_runs)
}

generate_permuted_rates <- function(genotype,genotype_asr_objs,permutations,df,tr,tip_names){
  # Model data
  comparitor_index_mat <- genotype_asr_objs[[genotype]]$index.mat
  comparitor_p <- sapply(1:max(comparitor_index_mat, na.rm = TRUE), function(x)
    na.omit(c(genotype_asr_objs[[genotype]]$solution))[na.omit(c(comparitor_index_mat) == x)][1])
  
  reconstructed_parent_child_df <- lapply(permutations[[genotype]],FUN=function(x){
    dataset <- cbind(tip_names,x)
    outcome_str <- setNames(x,tip_names)
    recon_result <- ancRECON(phy = tr,data = dataset,p = comparitor_p, method = 'joint', rate.cat = 1, rate.mat = comparitor_index_mat, root.p = comparitor_p, get.likelihood = F, get.tip.states = F)
    asr_result <- phyloAMR::get_parent_child_data(tr = tr, ancestral_states =  recon_result$lik.anc.states,trait_data =  outcome_str, confidence_threshold =  NULL,node_states = 'joint')
    asr_parent_child_df <- get_continuation_data(asr_result, node_states = 'joint')
  })
  
  return(reconstructed_parent_child_df)
}

generate_synchronous_test_pvalue <- function(trait_synchronous_obs,trait_parent_child,permutation_parent_child,num_permutations){
  comparitor_sychronous_perm <- lapply(permutation_parent_child,FUN=function(x){
    synchronous_transitions(x, trait_parent_child)
  })
  asr_permutation <- do.call(rbind,comparitor_sychronous_perm)
  
  # Calculate p-values for each event of interest and add them to the trait observation data
  transition_pval <- c(1 + sum(asr_permutation$synchronous_transitions_num >= trait_synchronous_obs$synchronous_transitions_num)) / c(1 + num_permutations)
  synchronous_gain_pval <- c(1 + sum(asr_permutation$synchronous_gains_num >= trait_synchronous_obs$synchronous_gains_num)) / c(1 + num_permutations)
  synchronous_loss_pval <- c(1 + sum(asr_permutation$synchronous_losses_num >= trait_synchronous_obs$synchronous_losses_num)) / c(1 + num_permutations)
  synchronous_gain_loss_pval <- c(1 + sum(asr_permutation$synchronous_gain_loss_num >= trait_synchronous_obs$synchronous_gain_loss_num)) / c(1 + num_permutations)
  synchronous_loss_gain_pval <- c(1 + sum(asr_permutation$synchronous_loss_gain_num >= trait_synchronous_obs$synchronous_loss_gain_num)) / c(1 + num_permutations)
  
  results_df <- data.frame(transition_pval,synchronous_gain_pval,synchronous_loss_pval,synchronous_gain_loss_pval,synchronous_loss_gain_pval)
  
  results <- list(results_df = results_df,syncronous_permutation = asr_permutation)
  return(results)
}

generate_downstream_test_pvalue <- function(trait_downstream_result,trait_parent_child,permutation_parent_child,stretches,num_permutations){
  comparitor_downstream_perm <- lapply(permutation_parent_child,FUN=function(x){
    phyloAMR:::get_gain_loss_on_stretches(comparitor_parent_child_df = x, downstream_nodes =  stretches, node_states = 'joint', confidence = NULL)
  })
  asr_permutation <- do.call(rbind,comparitor_downstream_perm)
  
  # Calculate p-values for each event of interest and add them to the trait observation data
  downstream_transition_pval <-  c(1 + sum(asr_permutation$transitions_num >= trait_downstream_result$transitions_num)) / c(1 + num_permutations)
  downstream_transition_stretches_pval <-  c(1 + sum(asr_permutation$stretches_w_transitions_num >= trait_downstream_result$stretches_w_transitions_num)) / c(1 + num_permutations)
  downstream_gain_pval <- c(1 + sum(asr_permutation$gains_num >= trait_downstream_result$gains_num)) / c(1 + num_permutations)
  downstream_gains_stretches_pval <- c(1 + sum(asr_permutation$stretches_w_gains_num >= trait_downstream_result$stretches_w_gains_num)) / c(1 + num_permutations)
  downstream_loss_pval <- c(1 + sum(asr_permutation$loss_num >= trait_downstream_result$loss_num)) / c(1 + num_permutations)
  downstream_loss_stretches_pval <- c(1 + sum(asr_permutation$stretches_w_losses_num >= trait_downstream_result$stretches_w_losses_num)) / c(1 + num_permutations)
  
  results_df <- data.frame(downstream_transition_pval,downstream_transition_stretches_pval,downstream_gain_pval,downstream_gains_stretches_pval,downstream_loss_pval,downstream_loss_stretches_pval)
  
  results <- list(results_df = results_df,downstream_permutations = asr_permutation)
  return(results)
}