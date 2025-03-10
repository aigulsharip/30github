# Load necessary libraries
install.packages("glmnet")  # Install if not already installed
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(MASS)
library(glmnet)


###############################################################################
### Data Preprocessing step
##############################################################################
# Load the dataset
df <- read_excel("gwas_AD_data_metadata_no_empty.xlsx", sheet = "filtered_id_data")

# Check dataset
summary(df)
str(df)
dim(df)
View(df)
# Convert target variable to binary (1 = AD, 0 = Control)
df$group <- ifelse(df$group == "AD", 1, 0)

# Drop ID column if it exists
if ("id" %in% names(df)) {
  df <- df[, !names(df) %in% "id"]
}

###########################################
#Manual Feature selection
# Remove Non-Kazakh Samples
df <- df %>% filter(race == "kazakh")  # Keep only Kazakh samples

###########################################
# Handle Missing Data
# Count occurrences of "NA" as a string in each column BEFORE numeric conversion
na_string_counts <- sapply(df, function(x) sum(x == "NA", na.rm = TRUE))
print("Count of 'NA' strings per column:")
print(na_string_counts)

# Define variables for downstream analysis
exclude_vars <- c("group", "gender", "region", "race", "bmi")   #bmi removal to avoid redundancy
numeric_features <- setdiff(names(df), exclude_vars)  # Numeric variables

# Identify numeric columns with ≤ 30% missing data
missing_percent <- (na_string_counts[numeric_features]) / nrow(df)
numeric_features <- numeric_features[missing_percent <= 0.30]  # Filter numeric vars

# Convert numeric columns to numeric type
df[numeric_features] <- lapply(df[numeric_features], as.numeric)

# Handle missing values using mean imputation (only for numeric variables)
for (col in numeric_features) {
  df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
}

# Convert categorical variables to factors
df$gender <- as.factor(df$gender)  # Keep 'gender' as categorical
categorical_features <- c("gender")  # Gender is categorical but included in analysis

#Final feature list for analysis
analysis_features <- c( categorical_features, numeric_features)  # Combined for modeling

# Final variable sets for downstream analysis
print(paste("Numeric Features:", paste(numeric_features, collapse = ", ")))
print(paste("Categorical Features:", paste(categorical_features, collapse = ", ")))

# Confirm changes
analysis_features

###############################################################################
### Comparing Means: t-test or Wilcoxon test
###############################################################################

test_results <- lapply(numeric_features, function(var) {
  if (all(is.na(df[[var]]))) return(NULL)
  ad_values <- df[[var]][df$group == 1]
  control_values <- df[[var]][df$group == 0]
  
  if (shapiro.test(ad_values)$p.value > 0.05 & shapiro.test(control_values)$p.value > 0.05) {
    test <- t.test(ad_values, control_values, var.equal = FALSE)
  } else {
    test <- wilcox.test(ad_values, control_values)
  }
  
  return(data.frame(Variable = var, p_value = test$p.value))
})

# Combine results
test_results_df <- do.call(rbind, test_results)
print(test_results_df)

# Filter significant variables (p-value < 0.05)
significant_vars <- test_results_df %>% filter(p_value < 0.05)
print("Significant variables:")
print(significant_vars)

# Compute means and standard deviations for significant variables and save as a data frame
if (nrow(significant_vars) > 0) {
  summary_stats <- df %>%
    group_by(group) %>%
    summarise(across(all_of(significant_vars$Variable), list(mean = ~mean(., na.rm = TRUE), 
                                                             sd = ~sd(., na.rm = TRUE))))
  
  # Save statistics as a dataset in the environment
  significant_stats <- summary_stats
  
  # Print the table
  print("\nTable of Mean and Standard Deviation for Significant Variables:")
  print(significant_stats)
}


###########################################
### Lasso Logistic Regression
###########################################

#Prepare Data for Lasso Logistic Regression
# Extract feature matrix X and response variable y for FULL feature set
X_full <- as.matrix(df[, analysis_features])  # Convert to matrix for glmnet
y_full <- df$group  # Response variable

# Split into training and test sets
set.seed(123)
train_index <- createDataPartition(y_full, p = 0.8, list = FALSE)
X_train_full <- X_full[train_index, ]
X_test_full <- X_full[-train_index, ]
y_train_full <- y_full[train_index]
y_test_full <- y_full[-train_index]

###########################################
### Perform Lasso Logistic Regression on Full Feature Set

# Train Lasso model (L1 regularization)
lasso_model_full <- glmnet(X_train_full, y_train_full, alpha = 1, family = "binomial")

# Cross-validation to find best lambda
set.seed(123)
cv_lasso_full <- cv.glmnet(X_train_full, y_train_full, alpha = 1, family = "binomial")

# Extract optimal lambda
best_lambda_full <- cv_lasso_full$lambda.min
print(paste("Optimal lambda (full feature set):", best_lambda_full))

# Fit final model with best lambda
lasso_final_full <- glmnet(X_train_full, y_train_full, alpha = 1, family = "binomial", lambda = best_lambda_full)

###########################################
### Identify Important Features (Full Feature Set)
# Extract feature coefficients
lasso_coef_full <- coef(lasso_final_full)

# Convert to dataframe
coef_df_full <- data.frame(
  Feature = rownames(lasso_coef_full),
  Coefficient = as.vector(lasso_coef_full)
) %>%
  filter(Feature != "(Intercept)", Coefficient != 0) %>%
  arrange(desc(abs(Coefficient)))

# Print selected features
print("Selected predictive variables (full feature set):")
print(coef_df_full)

###########################################
### Visualize Feature Importance (Full Feature Set)

ggplot(coef_df_full, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Feature Importance (Lasso Logistic Regression - Full Feature Set)",
       x = "Biomarker",
       y = "Coefficient",
       fill = "Effect Direction")

###########################################
### Model Evaluation (Full Feature Set)

# Make predictions
pred_probs_full <- predict(lasso_final_full, newx = X_test_full, type = "response")
pred_classes_full <- ifelse(pred_probs_full > 0.5, 1, 0)

# Create confusion matrix
conf_matrix_full <- confusionMatrix(factor(pred_classes_full), factor(y_test_full))

# Extract precision, recall, and F1-score
precision_full <- conf_matrix_full$byClass["Precision"]
recall_full <- conf_matrix_full$byClass["Sensitivity"]
f1_score_full <- 2 * ((precision_full * recall_full) / (precision_full + recall_full))

# Print results
print(paste("Precision (full feature set):", round(precision_full, 3)))
print(paste("Recall (Sensitivity, full feature set):", round(recall_full, 3)))
print(paste("F1-score (full feature set):", round(f1_score_full, 3)))

# Compute and plot ROC curve
roc_curve_full <- roc(y_test_full, pred_probs_full)
auc_value_full <- auc(roc_curve_full)

# Print AUC-ROC
print(paste("AUC-ROC (full feature set):", round(auc_value_full, 3)))

# Plot ROC Curve
plot(roc_curve_full, col="blue", lwd=2, main="ROC Curve (Full Feature Set)")
abline(a=0, b=1, lty=2, col="gray")  # Add diagonal line


###########################################
### Lasso Logistic Regression without Cognitive Scores (clockt, mmse)
###########################################
# Remove cognitive variables
df_no_cognitive <- df %>% select(-clockt, -mmse)

# Extract feature matrix X and response variable y
X_no_cognitive <- as.matrix(df_no_cognitive[, analysis_features[analysis_features %in% colnames(df_no_cognitive)]])  
y_no_cognitive <- df_no_cognitive$group  

# Split into training and test sets
set.seed(123)
train_index <- createDataPartition(y_no_cognitive, p = 0.8, list = FALSE)
X_train_no_cognitive <- X_no_cognitive[train_index, ]
X_test_no_cognitive <- X_no_cognitive[-train_index, ]
y_train_no_cognitive <- y_no_cognitive[train_index]
y_test_no_cognitive <- y_no_cognitive[-train_index]

# Perform Lasso Logistic Regression
lasso_model_no_cognitive <- glmnet(X_train_no_cognitive, y_train_no_cognitive, alpha = 1, family = "binomial")

# Cross-validation to find best lambda
set.seed(123)
cv_lasso_no_cognitive <- cv.glmnet(X_train_no_cognitive, y_train_no_cognitive, alpha = 1, family = "binomial")
best_lambda_no_cognitive <- cv_lasso_no_cognitive$lambda.min
print(paste("Optimal lambda (without cognitive scores):", best_lambda_no_cognitive))

# Fit final model with best lambda
lasso_final_no_cognitive <- glmnet(X_train_no_cognitive, y_train_no_cognitive, alpha = 1, family = "binomial", lambda = best_lambda_no_cognitive)

###########################################
### Identify Important Features (Without Cognitive Scores)

# Extract feature coefficients
lasso_coef_no_cognitive <- coef(lasso_final_no_cognitive)

# Convert to dataframe
coef_df_no_cognitive <- data.frame(
  Feature = rownames(lasso_coef_no_cognitive),
  Coefficient = as.vector(lasso_coef_no_cognitive)
) %>%
  filter(Feature != "(Intercept)", Coefficient != 0) %>%
  arrange(desc(abs(Coefficient)))

# Print selected features
print("Selected predictive variables (without cognitive scores):")
print(coef_df_no_cognitive)

###########################################
### Visualize Feature Importance (Without Cognitive Scores)

ggplot(coef_df_no_cognitive, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Feature Importance Without Cognitive Scores (Lasso Regression)",
       x = "Biomarker",
       y = "Coefficient",
       fill = "Effect Direction")

###########################################
### Model Evaluation (Without Cognitive Scores)

# Make predictions
pred_probs_no_cognitive <- predict(lasso_final_no_cognitive, newx = X_test_no_cognitive, type = "response")
pred_classes_no_cognitive <- ifelse(pred_probs_no_cognitive > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_no_cognitive <- confusionMatrix(factor(pred_classes_no_cognitive), factor(y_test_no_cognitive))
precision_no_cognitive <- conf_matrix_no_cognitive$byClass["Precision"]
recall_no_cognitive <- conf_matrix_no_cognitive$byClass["Sensitivity"]
f1_score_no_cognitive <- 2 * ((precision_no_cognitive * recall_no_cognitive) / (precision_no_cognitive + recall_no_cognitive))

print(paste("Precision (without cognitive scores):", round(precision_no_cognitive, 3)))
print(paste("Recall (Sensitivity, without cognitive scores):", round(recall_no_cognitive, 3)))
print(paste("F1-score (without cognitive scores):", round(f1_score_no_cognitive, 3)))

# AUC-ROC
roc_curve_no_cognitive <- roc(y_test_no_cognitive, as.vector(pred_probs_no_cognitive))
auc_value_no_cognitive <- auc(roc_curve_no_cognitive)
print(paste("AUC-ROC (without cognitive scores):", round(auc_value_no_cognitive, 3)))

# Plot ROC Curve
plot(roc_curve_no_cognitive, col="blue", lwd=2, main="ROC Curve Without Cognitive Scores")
abline(a=0, b=1, lty=2, col="gray")  

###############################################################################
### Random Forest for predicting AD status
##############################################################################

# Install necessary packages if not installed
install.packages("randomForest")
install.packages("pROC")  # For AUC-ROC
install.packages("caret") # For confusion matrix
library(randomForest)
library(pROC)
library(caret)

# Prepare Data
X_full <- as.matrix(df[, numeric_columns])  # Use all numeric predictors
y_full <- df$group  # Response variable

# Split into training and test sets
set.seed(123)
train_index <- createDataPartition(y_full, p = 0.8, list = FALSE)
X_train_full <- X_full[train_index, ]
X_test_full <- X_full[-train_index, ]
y_train_full <- y_full[train_index]
y_test_full <- y_full[-train_index]

# Convert to data frame for randomForest
train_data <- data.frame(X_train_full, group = as.factor(y_train_full))
test_data <- data.frame(X_test_full, group = as.factor(y_test_full))

# Train Random Forest Model
set.seed(123)
rf_model <- randomForest(group ~ ., data = train_data, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model)

# Feature Importance
importance_df <- data.frame(Feature = rownames(rf_model$importance), Importance = rf_model$importance[, 1])
importance_df <- importance_df %>% arrange(desc(Importance))

# Print top important features
print("Top Predictive Features (Random Forest):")
print(importance_df[1:20, ])  # Top 20 features

# Visualize Feature Importance
ggplot(importance_df[1:20, ], aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col(fill = "blue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Feature Importance (Random Forest)",
       x = "Feature",
       y = "Importance Score")

# Make Predictions
rf_predictions <- predict(rf_model, newdata = test_data, type = "prob")[, 2]  # Get probability of AD
rf_pred_classes <- ifelse(rf_predictions > 0.5, 1, 0)  # Convert to binary

# Confusion Matrix
conf_matrix_rf <- confusionMatrix(factor(rf_pred_classes), factor(y_test_full))
precision_rf <- conf_matrix_rf$byClass["Precision"]
recall_rf <- conf_matrix_rf$byClass["Sensitivity"]
f1_score_rf <- 2 * ((precision_rf * recall_rf) / (precision_rf + recall_rf))

print(paste("Precision (Random Forest):", round(precision_rf, 3)))
print(paste("Recall (Sensitivity, Random Forest):", round(recall_rf, 3)))
print(paste("F1-score (Random Forest):", round(f1_score_rf, 3)))

# AUC-ROC
roc_curve_rf <- roc(y_test_full, rf_predictions)
auc_value_rf <- auc(roc_curve_rf)
print(paste("AUC-ROC (Random Forest):", round(auc_value_rf, 3)))

# Plot ROC Curve
plot(roc_curve_rf, col="blue", lwd=2, main="ROC Curve (Random Forest)")
abline(a=0, b=1, lty=2, col="gray")  # Add diagonal line

