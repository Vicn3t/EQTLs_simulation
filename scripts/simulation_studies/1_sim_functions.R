#######################################################################################################
#                                                                                                     #
#                                    SIMULATION DATA PREPARATTION                                     #
#                                                                                                     #
#######################################################################################################

compute_variants_maf <- function(genotypes, n_variants_to_sample, breeds_prefix) {
  
  # Extract positions from rownames
  variant_ids <- rownames(genotypes)
  variant_order <- str_extract(variant_ids, "(?<=_)\\d+$")
  variant_position <- str_extract(variant_ids, "(?<=_)\\d+(?=_)")
  
  # Remove rows with NA from extraction
  valid_rows <- !is.na(variant_order) & !is.na(variant_position)
  genotypes <- genotypes[valid_rows, , drop = FALSE]
  variant_ids <- variant_ids[valid_rows]
  variant_order <- variant_order[valid_rows]
  variant_position <- variant_position[valid_rows]
  
  # Compute MAF
  variants_metadata <- data.frame(
    variant_order    = as.integer(variant_order),
    variant_id       = variant_ids,
    variant_position = as.integer(variant_position),
    maf              = rowSums(genotypes) / (2 * ncol(genotypes))
  )
  
  # Compute breed-specific MAFs
  for (breed in breeds_prefix) {
    breed_inds <- grepl(breed, colnames(genotypes))
    variants_metadata[[paste0(breed, "_maf")]] <-
      rowSums(genotypes[, breed_inds, drop = FALSE]) / (2 * sum(breed_inds))
  }
  
  # Filter out extreme MAFs
  variants_metadata <- variants_metadata %>%
    filter(maf > 0.05 & maf < 0.95)

  
  # Add MAF-based category
  variants_metadata <- variants_metadata %>%
    mutate(variant_category = case_when(
      !!sym(paste0(breeds_prefix[1], "_maf")) == 0 ~ paste0(breeds_prefix[2], "_specific"),
      !!sym(paste0(breeds_prefix[2], "_maf")) == 0 ~ paste0(breeds_prefix[1], "_specific"),
      (!!sym(paste0(breeds_prefix[1], "_maf")) - !!sym(paste0(breeds_prefix[2], "_maf"))) > 0.20 ~ paste0(breeds_prefix[1], "_contrasted"),
      (!!sym(paste0(breeds_prefix[2], "_maf")) - !!sym(paste0(breeds_prefix[1], "_maf"))) > 0.20 ~ paste0(breeds_prefix[2], "_contrasted"),
      TRUE ~ "homogeneous"
    ))
  
  # Init columns
  variants_metadata$breed_effect <- "null_effect"
  variants_metadata$sim_beta     <- 0
  for (breed in breeds_prefix) {
    variants_metadata[[paste0(breed, "_sim_beta")]] <- 0
  }
  
  return(variants_metadata)
  
}

#######################################################################################################

### Generate multiple plots related to MAF distributions
MAF_plots <- function(variants_metadata, breeds_prefix) {

  # Create output directory if it doesn't exist
  output_dir <- here("data/04_simulation_r_data/MAF_plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- Histogram: Global MAF distribution ---
  p1 <- ggplot(variants_metadata, aes(x = maf)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, fill = "steelblue", 
                   color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Global MAF Distribution", x = "Global MAF", y = "Density")
  
  ggsave(file.path(output_dir, "maf_Distribution.png"), plot = p1, width = 8, height = 6, dpi = 300)
  plot(p1)
  
  # --- Barplot: Variant distribution by category ---
  p2 <- ggplot(variants_metadata, aes(x = variant_category, fill = variant_category)) +
    geom_bar() +
    theme_minimal() +
    labs(title = "Variant Distribution by Category", x = "Category", y = "Number of Variants") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = "none")
  
  ggsave(file.path(output_dir, "variants_Distribution_by_Category.png"), plot = p2, 
         width = 8, height = 6, dpi = 300)
  plot(p2)
  
  # --- Density plot: MAF by breed (dynamic) ---
  # Gather breed-specific MAFs into long format
  maf_long <- variants_metadata %>%
    select(all_of(paste0(breeds_prefix, "_maf"))) %>%
    rename_with(~ str_remove(., "_maf"), everything()) %>%
    mutate(variant_id = rownames(variants_metadata)) %>%
    pivot_longer(
      cols = all_of(breeds_prefix),
      names_to = "Breed", values_to = "MAF")
  
  p3 <- ggplot(maf_long, aes(x = MAF, fill = Breed)) +
    geom_density(alpha = 0.5, color = "black") +
    theme_minimal() +
    labs(title = "MAF Distribution by Breed", x = "MAF", y = "Density") +
    guides(fill = guide_legend(title = "Breed"))
  
  ggsave(file.path(output_dir, "maf_Distribution_by_Breed.png"), plot = p3, 
         width = 8, height = 6, dpi = 300)
  plot(p3)

}


#######################################################################################################

### Assign n variants per gene within ±10Mb of the gene center
associate_genes_with_qtls <- function(genes_metadata, variants_metadata, 
                                      variants_gene_max_distance, n_qtls_per_gene) {
  
  # For each gene, assign closest variant positions within the defined window
  genes_metadata$associated_variant_pos <- lapply(1:nrow(genes_metadata), function(i) {
    gene <- genes_metadata[i, ]
    
    # Compute gene center and QTL window
    gene_center  <- (gene$start + gene$end) / 2
    lower_bound  <- gene_center - variants_gene_max_distance
    upper_bound  <- gene_center + variants_gene_max_distance
    
    # Filter variants within the window
    surrounding_variants <- variants_metadata %>%
      filter(variant_position > lower_bound, variant_position < upper_bound)
    
    # Sample n QTLs per gene or return NA if not enough
    if (nrow(surrounding_variants) >= n_qtls_per_gene) {
      sample(surrounding_variants$variant_position, n_qtls_per_gene)
    } else {
      rep(NA, n_qtls_per_gene)
    }
  })
  
  # Build a QTL-gene association table (one row = one QTL-gene link)
  gene_qtl_association <- data.frame(
    qtl_position    = unlist(genes_metadata$associated_variant_pos),
    associated_gene = rep(genes_metadata$gene_id, 
                          times = sapply(genes_metadata$associated_variant_pos, length))
  )
  
  # Find corresponding rows in variants_metadata
  variant_table_qtl_id <- match(gene_qtl_association$qtl_position, variants_metadata$variant_position)
  qtl_table <- variants_metadata[variant_table_qtl_id, ]
  
  # Combine gene ID with associated QTL metadata
  gene_qtl_associations_table <- cbind(
    gene_qtl_association[, "associated_gene", drop = FALSE],
    qtl_table
  )
  
  gene_qtl_associations_table <- gene_qtl_associations_table %>%
    arrange(associated_gene, variant_order)
  return(gene_qtl_associations_table)
}

#######################################################################################################

simulate_qtl_effects <- function(gene_qtl_associations_table, breeds_prefix, n_qtls_per_gene, 
                                 breed_effects_proportions, mean_effect, var_effect) {
  
  # Working directly on the gene qtl association table
  working_qtls_df = gene_qtl_associations_table
  
  qtl_number = nrow(working_qtls_df)
  
  # Assign a global beta from normal distribution
  working_qtls_df$sim_beta = sample(c(1, -1), qtl_number, replace = TRUE) * 
    rnorm(qtl_number, mean_effect, var_effect)
  
  
  # Find potential specific QTLs of both breeds
  for (breed in breeds_prefix) {
    
    # Fetch the breed beta column name and check for breed specific QTLs
    breed_col <- paste0(breed, "_sim_beta")
    specific_idx <- which(grepl(paste0(breed, "_specific"), working_qtls_df$variant_category))
    
    if (length(specific_idx) > 0) {
      working_qtls_df[[breed_col]][specific_idx] <- 
        sample(c(1, -1), length(specific_idx), replace = TRUE) * 
        rnorm(length(specific_idx), mean_effect, var_effect)
      working_qtls_df$breed_effect[specific_idx] <- "specific_effect"
    }
  }
  
  # Identify indices of QTLs that are not breed-specific
  idx <- which(working_qtls_df$breed_effect == "null_effect")
  
  
  # Define in the function call the proportions for the three effect types (used only for non-specific QTLs) :
  # - same_effect: same effect in both breeds
  # - different_effect: different effect sizes between breeds
  # - opposite_effect: effects in opposite directions between breeds
  
  # Randomly assign an effect type to each non-specific QTL according to the defined proportions
  working_qtls_df$breed_effect[idx] <- sample(names(breed_effects_proportions), length(idx), 
                                              replace = TRUE, prob = breed_effects_proportions)
  
  # Get the effect types for those QTLs
  eff_types <- working_qtls_df$breed_effect[idx]
  
  # Prepare empty vectors for simulated betas
  n <- length(idx)
  beta1 <- numeric(n)
  beta2 <- numeric(n)
  
  # Identify rows by effect type
  same_idx <- eff_types == "same_effect"
  diff_idx <- eff_types == "different_effect"
  opp_idx  <- eff_types == "opposite_effect"
  
  # Generate beta1 and beta2 for each effect type
  beta1[same_idx] <- rnorm(sum(same_idx), mean_effect, var_effect) * sample(c(-1, 1), sum(same_idx), replace = TRUE)
  beta2[same_idx] <- beta1[same_idx]
  
  common_signs <- sample(c(-1, 1), sum(diff_idx), replace = TRUE)
  beta1[diff_idx] <- rnorm(sum(diff_idx), mean_effect * 0.5, var_effect) * common_signs
  beta2[diff_idx] <- rnorm(sum(diff_idx), mean_effect * 1.5, var_effect) * common_signs
  
  beta1[opp_idx] <- rnorm(sum(diff_idx), mean_effect * 1, var_effect)
  beta2[opp_idx] <- rnorm(sum(diff_idx), mean_effect * -1, var_effect)
  
  # Swap half of the beta values between the two breeds to avoid breed-specific bias
  # (e.g. : avoid breed1 systematically getting higher/lower effects)
  swap_indices <- sample(n, size = floor(n / 2))
  
  # Perform the swap
  temp <- beta1[swap_indices]
  beta1[swap_indices] <- beta2[swap_indices]
  beta2[swap_indices] <- temp
  
  
  breed1 <- breeds_prefix[1]
  breed2 <- breeds_prefix[2]
  
  # Store the simulated beta values in the dataframe for both breeds
  working_qtls_df[[paste0(breed1, "_sim_beta")]][idx] <- beta1
  working_qtls_df[[paste0(breed2, "_sim_beta")]][idx] <- beta2
  
  saveRDS(working_qtls_df, here("data/04_simulation_r_data/gene_qtl_association_beta_table.RDS"))
  
  return(working_qtls_df)
}

#######################################################################################################


simulate_gene_expression <- function(qtls_beta_table, genotypes, breeds_genotypes, breeds_prefix, 
                                     variants_metadata, heritability) {
  
  # Loop through each unique gene
  gene_expression_list <- lapply(unique(qtls_beta_table$associated_gene), function(gene_id) {
    
    # Filter QTL info for the current gene
    gene_qtls_df <- subset(qtls_beta_table, associated_gene == gene_id)
    
    # Create a working copy of the variant MAF table
    variant_table <- variants_metadata
    
    # Update simulated beta values
    idx <- match(gene_qtls_df$variant_id, variant_table$variant_id)
    variant_table[idx,] <- gene_qtls_df[, names(gene_qtls_df) != "associated_gene"]
    
    # Global genetic value computation
    beta_matrix <- matrix(variant_table$sim_beta, ncol = 1, dimnames = list(variant_table$variant_id, "sim_beta"))
    genetic_value <- as.vector(t(genotypes) %*% beta_matrix)
    
    # Add environmental noise while controlling heritability
    error_var <- (var(genetic_value) * (1 - heritability)) / heritability
    noise <- rnorm(length(genetic_value), mean = 0, sd = sqrt(error_var))
    expression <- genetic_value + noise
    
    # Check resulting heritability
    h2_check <- var(genetic_value) / (var(genetic_value) + var(noise))
    
    # Format global expression output
    global_expr_df <- data.frame(
      id = colnames(genotypes),
      expr = expression
    )
    colnames(global_expr_df)[2] <- paste0(gene_id, "_global_expression")
    
    # Intra-breed expression simulation
    breed_gene_expression_list <- lapply(breeds_prefix, function(breed) {
      breed_beta_matrix <- matrix(
        variant_table[[paste0(breed, "_sim_beta")]],
        ncol = 1,
        dimnames = list(variant_table$variant_id, paste0(breed, "_sim_beta"))
      )
      
      breed_genetic_value <- as.vector(t(breeds_genotypes[[breed]]) %*% breed_beta_matrix)
      breed_error_var <- (var(breed_genetic_value) * (1 - heritability)) / heritability
      breed_noise <- rnorm(length(breed_genetic_value), mean = 0, sd = sqrt(breed_error_var))
      breed_expression <- breed_genetic_value + breed_noise
      
      h2_breed <- var(breed_genetic_value) / (var(breed_genetic_value) + var(breed_noise))
      
      # Format breed expression output
      breed_expr_df <- data.frame(
        id = colnames(breeds_genotypes[[breed]]),
        expr = breed_expression,
        stringsAsFactors = FALSE
      )
      colnames(breed_expr_df) <- c("id", paste0(gene_id, "_breeds_expression"))
      
      return(breed_expr_df)
    })
    
    # Combine the breed-specific simulated expression data frames
    breed_expression_df <- dplyr::bind_rows(breed_gene_expression_list)
    
    merged_df <- merge(global_expr_df, breed_expression_df, by = "id")
    
    return(merged_df)
    
  })
  
  gene_expression_list <- setNames(gene_expression_list,unique(qtls_beta_table$associated_gene) )
  
  return(gene_expression_list)
}


#######################################################################################################

# Function to extract global expression matrix from gene_expression_lists
get_global_gene_expression_df <- function(gene_expression_lists) {
  
  # Extract the global expression data (2nd column of each list's first element)
  global_expr_matrix <- lapply(gene_expression_lists, function(gene_entry) {
    expr_values <- gene_entry[, 2, drop = FALSE]  # 2nd column
    rownames(expr_values) <- gene_entry[, 1] # use IDs as row names
    return(expr_values)
  })
  
  # Combine all into a single data frame: columns = genes, rows = individuals
  global_gene_expression_df <- do.call(cbind, global_expr_matrix)
  colnames(global_gene_expression_df) <- names(gene_expression_lists)
  
  saveRDS(global_gene_expression_df, here("data/04_simulation_r_data","global_gene_expression_df.RDS"))
  
  return(global_gene_expression_df)
}

# Function to extract breed-specific expression matrix from gene_expression_lists
# Returns a named list of data frames: one per breed
get_breeds_gene_expression_df <- function(gene_expression_lists) {
  
  # Extract the breeds expression data (2nd column of each list's first element)
  breeds_expr_matrix <- lapply(gene_expression_lists, function(gene_entry) {
    expr_values <- gene_entry[, 3, drop = FALSE]  # 3nd column
    rownames(expr_values) <- gene_entry[, 1] # use IDs as row names
    return(expr_values)
  })
  
  # Combine all into a single data frame: columns = genes, rows = individuals
  breeds_gene_expression_df <- do.call(cbind, breeds_expr_matrix)
  colnames(breeds_gene_expression_df) <- names(gene_expression_lists)
  
  saveRDS(breeds_gene_expression_df, here("data/04_simulation_r_data","breeds_gene_expression_df.RDS"))
  
  
  return(breeds_gene_expression_df)
}

#######################################################################################################
#                                                                                                     #
#                                              MATRIX_eQTL                                            #
#                                                                                                     #
#######################################################################################################


#' Create PCA covariates as SlicedData object
get_pca_covariates <- function(pca_rds_path, individual_ids) {
  
  # Load PCA result and extract PC coordinates (PCs x individuals)
  pcs <- t(readRDS(pca_rds_path)$ind$coord)
  
  # Reorder columns to match desired individual order
  pcs <- pcs[, individual_ids, drop = FALSE]
  
  # Create SlicedData object
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(pcs)
  
  return(cvrt)
}



#######################################################################################################

run_eqtl_model <- function(genotypes_matrix, expression_df, covariates_list, output_prefix, output_dir,Matrix_eQTL_verbose) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  sliced_variant <- SlicedData$new()
  sliced_variant$CreateFromMatrix(genotypes_matrix)
  
  alpha <- 0.05
  tests_number <- dim(expression_df)[2] * dim(genotypes_matrix)[1]  # gene_number * SNPs
  p_value_threshold <- alpha / tests_number
  
    sliced_expr <- SlicedData$new()
    sliced_expr$CreateFromMatrix(t(as.matrix(expression_df)))
    
    for (cvrt_name in names(covariates_list)) {
      cvrt_obj <- covariates_list[[cvrt_name]]
      model_name <- paste0(output_prefix, "_", cvrt_name)
      message("Running Matrix_eQTL: ", model_name)
      
      result <- Matrix_eQTL_main(
        snps = sliced_variant,
        gene = sliced_expr,
        cvrt = cvrt_obj,
        output_file_name = "",
        pvOutputThreshold = p_value_threshold,
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = Matrix_eQTL_verbose
      )$all$eqtls %>%
        mutate(variant_id = snps) %>%
        select(gene, variant_id, beta, pvalue)
      
      saveRDS(result, file = file.path(output_dir, paste0(model_name, ".RDS")))
    }
  }



#######################################################################################################


run_intra_breed_eqtl <- function(breed, genotypes, expression_data, metadata, add_no_cov_model, Matrix_eQTL_verbose) {
  
  ids <- metadata$id[metadata$prefix == breed]
  
  breed_genotypes <- genotypes[, ids]
  
  breed_genotypes <- breed_genotypes[sapply(as.data.frame(t(breed_genotypes)), var, na.rm = TRUE) > 0, ]
  
  breed_expression <- expression_data[ids, ]
  
  pca_path <- here("data/03_cleaned_r_data", paste0("geno_PCA_", breed, ".RDS"))
  breed_pca <- get_pca_covariates(pca_rds_path = pca_path, individual_ids = ids)
  
  if (add_no_cov_model){
    breed_covs_list <- list( pca_cov = breed_pca,no_cov = SlicedData$new())
  }else{breed_covs_list <- list(pca_cov = breed_pca)}
  
  run_eqtl_model(
    genotypes_matrix = breed_genotypes,
    expression_df = breed_expression,
    covariates_list = breed_covs_list,
    output_prefix = paste0("intra_", breed),
    output_dir = here("data", "05_sim_models_results", "Matrix_eQTL"),
    Matrix_eQTL_verbose = Matrix_eQTL_verbose
  )
}


#######################################################################################################
#                                                                                                     #
#                                   MODELS PERFORMANCE EVALUATION                                     #
#                                                                                                     #
#######################################################################################################

# Not used anymore

# bin_results <- function(df, calls_bin_size) {
#   df$bin <- NA_integer_
#   
#   for (gene in unique(df$gene_id)) {
#     sub_df <- df %>% filter(gene_id == gene)
#     bin_index <- 1
#     for (i in seq_len(nrow(sub_df))) {
#       if (i == 1 || (sub_df$variant_pos[i] - sub_df$variant_pos[i - 1]) > calls_bin_size) bin_index <- bin_index + (i != 1)
#       sub_df$bin[i] <- bin_index
#     }
#     df$bin[df$gene_id == gene] <- sub_df$bin
#   }
#   
#   df %>% group_by(gene_id, bin) %>%
#     ungroup()
# }

#######################################################################################################

### Results data prep: cleaning, merging and annotating model results vs simulated QTLs

# Step 1: Clean model results
# - Rename columns for consistency
# - Extract variant position from variant ID
prepare_model_results <- function(model_df) {
  model_df %>%
    rename(gene_id = gene, estimated_beta = beta, adjusted_pvalue = pvalue) %>%
    mutate(variant_pos = as.numeric(str_extract(variant_id, "(?<=_)\\d+(?=_)"))) %>%
    arrange(gene_id, variant_pos)
}

# Step 2: Optionally bin the variants
# - If enabled, select top variants per region (to reduce LD-driven redundancy)
apply_binning_if_needed <- function(df, binning = FALSE) {
  if (binning) bin_results(df, calls_bin_size = 5000) else df
}

# Step 3: Annotate detected calls with variant metadata
# - Join variant annotations
# - Mark variants as "detected QTLs"
annotate_detected_qtls <- function(df, all_variants_df) {
  df %>%
    left_join(all_variants_df, by = "variant_id") %>%
    mutate(is_detected_as_QTL = TRUE)
}

# Step 4: Merge model detections with simulated QTLs
# - Join simulated QTL data and annotated calls
# - Join is done on both QTL and variant characteristics to ensure accurate matching
merge_with_qtls <- function(detected_df, qtls_df, all_variants_df) {
  qtls_annotated <- qtls_df %>%
    left_join(all_variants_df, by = "variant_id")
  
  full_join(qtls_annotated, detected_df,
            by = c("gene_id", "variant_id", "variant_pos", "maf", "LD_maf", "LW_maf", "variant_category"))
}

# Step 5: Fill missing fields
# - For variants not matched between QTLs and calls, fill default values
# - Ensure consistent downstream logic (e.g., classification TP/FP/FN)
fill_missing_info <- function(df) {
  df %>%
    mutate(
      is_qtl = if_else(is.na(is_qtl), FALSE, is_qtl),
      is_detected_as_QTL = if_else(is.na(is_detected_as_QTL), FALSE, is_detected_as_QTL),
      breed_effect = if_else(is.na(breed_effect), "null_effect", breed_effect),
      sim_beta = if_else(is.na(sim_beta), 0, sim_beta),
      LD_sim_beta = if_else(is.na(LD_sim_beta), 0, LD_sim_beta),
      LW_sim_beta = if_else(is.na(LW_sim_beta), 0, LW_sim_beta),
      adjusted_pvalue = if_else(is.na(adjusted_pvalue), 0.05, adjusted_pvalue),
      estimated_beta = if_else(is.na(estimated_beta), 0, estimated_beta)
    )
}

# Step 6: Classify variants
# Classify each variant as TP, FP or FN based on QTL status and detection
add_detection_status <- function(df) {
  df %>%
    mutate(
      detection_status = case_when(
        is_qtl & is_detected_as_QTL  ~ "TP",
        !is_qtl & is_detected_as_QTL ~ "FP",
        is_qtl & !is_detected_as_QTL ~ "FN",
        TRUE                         ~ NA_character_
      )
    )
}

# Step 7: filter undetactable variants and eQTLs
# For intra models, filter variants with no variability (breed MAF = 0 or 1)
filter_no_variability_variants <- function(df, model_name, breeds_prefix) {
  
  for (prefix in breeds_prefix) {
    if (grepl(prefix, model_name)) {
      breed_maf_col <- paste0(prefix, "_maf")
      
      # Filter variants with MAF strictly between 0 and 1
      df <- df %>%
        filter(.data[[breed_maf_col]] > 0 & .data[[breed_maf_col]] < 1)
      
      break  # Stop loop once the relevant breed is found
    }
  }
  
  return(df)  # Return the filtered dataframe
}

# Step 6: Execute every function
# - Apply all steps in order to return a fully annotated merged table
merge_qtls_and_model_calls <- function(model_results_df, model_name, qtls_df, all_variants_df, is_binning_enabled,breeds_prefix) {
  clean_df    <- prepare_model_results(model_results_df)
  binned_df   <- apply_binning_if_needed(clean_df, is_binning_enabled)
  detected_df <- annotate_detected_qtls(binned_df, all_variants_df)
  merged_df   <- merge_with_qtls(detected_df, qtls_df, all_variants_df)
  filled_df    <- fill_missing_info(merged_df)
  detection_status_df    <- add_detection_status(filled_df)
  final_df <- filter_no_variability_variants(detection_status_df, model_name, breeds_prefix)
  
  return(final_df)
}

#######################################################################################################

# Process a single model result file: load data, merge with QTLs, and annotate with model name
process_model_file <- function(file_path, qtls_df, variants_df, is_binning_enabled, breeds_prefix) {
  
  # Load the model results (RDS file containing QTL detection results)
  model_results_df <- readRDS(file_path)
  
  # Extract model name from the file name (without extension)
  model_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Merge model calls with simulated QTLs and variant annotations
  result <- merge_qtls_and_model_calls(
    model_results_df = model_results_df,
    model_name = model_name,
    qtls_df = qtls_df,
    all_variants_df = variants_df,
    is_binning_enabled = is_binning_enabled,
    breeds_prefix = breeds_prefix
  )
  
  # Add model name to the result for downstream grouping/comparison
  result$model_name <- model_name
  
  result <- result %>%
    select(gene_id, variant_id, variant_pos, maf, LD_maf, LW_maf,variant_category, is_qtl, breed_effect, LD_sim_beta,
           LW_sim_beta, model_name, is_detected_as_QTL, detection_status, estimated_beta, adjusted_pvalue)
    
  return(result)
}

#######################################################################################################

# Browse all model result files and compute TP/FP/FN information for each model
compute_models_tp_fp_fn <- function(models_results_path, variants_df, qtls_df, is_binning_enabled = FALSE, breeds_prefix) {
  
  # Get all subdirectories (each corresponding to a different tool, e.g. Matrix_eQTL)
  subdirs <- list.dirs(models_results_path, recursive = FALSE)
  
  # Loop over subdirectories
  all_results <- lapply(subdirs, function(dir) {
    
    # Get all .RDS files (each one corresponding to a different model)
    files_paths <- list.files(dir, pattern = "\\.RDS$", full.names = TRUE)
    
    # Process each model result file
    lapply(files_paths, function(file_path) {
      process_model_file(file_path, qtls_df, variants_df, is_binning_enabled, breeds_prefix)
    })
  })
  
  # Combine all model results into a single data frame
  bind_rows(unlist(all_results, recursive = FALSE))
}



#######################################################################################################

# Compute performance metrics (TP, FP, FN, F1 Score, MCC) for each gene across models
evaluate_models_performances <- function(models_results_df, qtls_df, all_variants_df, breeds_prefix) {
  
  # Extract gene names (common across all models)
  genes_names <- unique(models_results_df$gene_id)
  
  # Extract model names
  models_names <- unique(models_results_df$model_name)
  
  # Extract all variant IDs (may be filtered by breed-specific MAF)
  variants_ids <- all_variants_df$variant_id
  
  # Loop over each model
  models_metrics_list <- lapply(models_names, function(model) {
    
    # Subset results for the current model
    model_results_df <- models_results_df %>%
      dplyr::filter(model_name == model)
    
    # Get QTLs for the current model
    qtls_metadata <- model_results_df %>%
      dplyr::filter(is_qtl == TRUE)
    
    qtls_ids <- qtls_metadata %>% 
      dplyr::pull(variant_id)
    
    # Filter variants based on breed-specific MAF if applicable
    for (prefix in breeds_prefix) {
      if (grepl(prefix, model)) {
        breed_maf_col <- paste0(prefix, "_maf")
        
        breed_variants_df <- all_variants_df %>%
          dplyr::filter(.data[[breed_maf_col]] > 0 & .data[[breed_maf_col]] < 1)
        
        variants_ids <- breed_variants_df %>%
          dplyr::pull(variant_id)
        break
      }
    }
    
    # Compute metrics for each gene
    metrics_list <- lapply(genes_names, function(gene) {
      
      model_calls <- model_results_df %>%
        dplyr::filter(gene_id == gene)
      
      TP <- sum(model_calls$detection_status == "TP", na.rm = TRUE)
      FP <- sum(model_calls$detection_status == "FP", na.rm = TRUE)
      FN <- sum(model_calls$detection_status == "FN", na.rm = TRUE)
      TN <- length(variants_ids) - TP - FP - FN
      
      # Compute precision and recall
      Precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      
      # Compute F1-score (set to 0 if both Precision and Recall are 0)
      F1_Score <-if ((Precision + Recall) > 0) {
        2 * (Precision * Recall) / (Precision + Recall)
        } else {
          0
        }
      
      # Compute Matthews Correlation Coefficient (MCC)
      denominator <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
      
      MCC <- if (denominator != 0) {
        ((TP * TN) - (FP * FN)) / denominator
      } else {
        0
      }
      
      # Return data frame with metrics for the current gene
      data.frame(
        model_name = model,
        gene_id    = gene,
        TP, FP, TN, FN,
        Recall, Precision, F1_Score, MCC
      )
    })
    
    # Combine gene-level metrics for the current model
    dplyr::bind_rows(metrics_list)
  })
  
  # Combine metrics for all models
  results<- dplyr::bind_rows(models_metrics_list)
  
  results <- results %>%
    dplyr::filter(!(TP == 0 & FN == 0)) 
  return(results)
  }


#######################################################################################################

# Summarize the model metrics
summarize_model_performances <- function(df) {
  
  results <- df %>%
    group_by(model_name) %>%
    summarise(
      n = n(),  # number of genes
      
      # Detection status
      TP = mean(TP),
      FP = mean(FP),
      TN = mean(TN),
      FN = mean(FN),
      
      # Precision + IC
      Precision_mean = mean(Precision),
      Precision_sd = sd(Precision),
      Precision_ic = 1.96 * Precision_sd / sqrt(n),
      
      # Recall + IC
      Recall_mean = mean(Recall),
      Recall_sd = sd(Recall),
      Recall_ic = 1.96 * Recall_sd / sqrt(n),
      
      # F1-score + IC
      F1_Score_mean = mean(F1_Score),
      F1_Score_sd = sd(F1_Score),
      F1_Score_ic = 1.96 * F1_Score_sd / sqrt(n),
      
      # MCC + IC
      MCC_mean = mean(MCC, na.rm = TRUE),
      MCC_sd = sd(MCC, na.rm = TRUE),
      MCC_ic = 1.96 * MCC_sd / sqrt(n),
      
      .groups = "drop"
    ) %>%
    arrange(desc(Recall_mean)) %>%
    select(-n)
  
  return(results)
}

#######################################################################################################

# we make the assumption that the intra simulation variability is low compared to the inter simulations
# variability
aggregate_simulation_metrics <- function(df) {
  df %>%
    group_by(model_name) %>%
    summarise(
      n = n(),
      
      # Detection status
      TP_mean = mean(TP),
      TP_sd = sd(TP),
      
      FP_mean = mean(FP),
      FP_sd = sd(TP),
      
      TN_mean = mean(TN),
      TN_sd = sd(TN),
      
      FN_mean = mean(FN),
      FN_sd = sd(FN),
      
      # Precision
      Precision = mean(Precision_mean, na.rm = TRUE),
      Precision_sd = sd(Precision_mean, na.rm = TRUE),
      Precision_ic = 1.96 * Precision_sd / sqrt(n),
      
      # Recall
      Recall = mean(Recall_mean, na.rm = TRUE),
      Recall_sd = sd(Recall_mean, na.rm = TRUE),
      Recall_ic = 1.96 * Recall_sd / sqrt(n),
      
      # F1-score
      F1_Score = mean(F1_Score_mean, na.rm = TRUE),
      F1_Score_sd = sd(F1_Score_mean, na.rm = TRUE),
      F1_Score_ic = 1.96 * F1_Score_sd / sqrt(n),
      
      # MCC
      MCC = mean(MCC_mean, na.rm = TRUE),
      MCC_sd = sd(MCC_mean, na.rm = TRUE),
      MCC_ic = 1.96 * MCC_sd / sqrt(n),
      
      .groups = "drop"
    )%>%
    select(-n)
}

#######################################################################################################

summarize_model_performances_long <- function(df) {
  df <- df %>%
    select(model_name, Precision, Recall, F1_Score,"MCC") %>%
    pivot_longer(
      cols = -model_name,
      names_to = "metric",
      values_to = "mean_value"
    )
  
  df$metric <- factor(df$metric,
                      levels = c("Recall","Precision", "F1_Score","MCC"))
  df
}

#######################################################################################################

