# Load libraries
library(glmnet)
library(glmnetUtils)
library(Metrics)
library(doParallel)
library(parallel)
library(foreach)

# Load training and testing datasets
load("training_filt")
load("testing_filt")

# Set number of cores based on computing capabilities
num_cores <- detectCores() - 32
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define loci and alpha values
loci_names <- setdiff(colnames(training_filt), "Age")
alpha_values <- c(0.000, 0.001, 0.008, 0.027, 0.064, 0.125, 0.216, 0.343, 0.512, 0.729, 1.000)

# Parallel loop using foreach with a specified number of iterations
results_list <- foreach(run = 1:600000, .packages = c("glmnet", "Metrics")) %dopar% {
  
  # 🔸 Log progress to per-worker file
  log_file <- paste0("progress_worker_", Sys.getpid(), ".log")
  cat(paste0("Run ", run, " started at ", Sys.time(), "\n"), 
      file = log_file, append = TRUE)
  
  # Sample 125 random loci
  selected_loci <- sample(loci_names, 125, replace = FALSE)
  
  # Prepare training and test matrices
  x_train <- as.matrix(training_filt[, selected_loci])
  y_train <- as.matrix(training_filt$Age)
  x_test  <- as.matrix(testing_filt[, selected_loci])
  y_test  <- as.matrix(testing_filt$Age)
  
  best_mae <- Inf
  best_alpha <- NA
  best_metrics <- NULL
  
  # Loop through alpha values
  for (a in alpha_values) {
    CVGLM <- cv.glmnet(
      x = x_train,
      y = y_train,
      nfolds = nrow(training_filt),
      alpha = a,
      type.measure = "mae",
      family = "gaussian",
      grouped = FALSE
    )
    
    # Predict on test data
    pred_test <- predict(CVGLM, newx = x_test, s = "lambda.min")
    pred_test <- exp(pred_test) - 1
    chronological_test <- exp(y_test) - 1
    mae_test_val <- mae(chronological_test, pred_test)
    mdae_test_val <- mdae(chronological_test, pred_test)
    
    # Track best alpha
    if (mae_test_val < best_mae) {
      best_mae <- mae_test_val
      best_alpha <- a
      
      pred_train <- predict(CVGLM, newx = x_train, s = "lambda.min")
      pred_train <- exp(pred_train) - 1
      chronological_train <- exp(y_train) - 1
      mae_train_val <- mae(chronological_train, pred_train)
      mdae_train_val <- mdae(chronological_train, pred_train)
      
      best_metrics <- data.frame(
        Run = run,
        Alpha = best_alpha,
        MAE_Train = mae_train_val,
        MDAE_Train = mdae_train_val,
        MAE_Test = mae_test_val,
        MDAE_Test = mdae_test_val,
        Loci = paste(selected_loci, collapse = ", ")
      )
    }
  }
  
  best_metrics
}

# Track progress
tail -f progress_worker_*.log

# Stop the cluster
stopCluster(cl)

# Combine into single data frame
results <- do.call(rbind, results_list)

# Save output
save(results, file="results")