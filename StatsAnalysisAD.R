# Load necessary libraries
library(readxl)
library(dplyr)
library(caret)

# Load Excel file
file_path <- "gwas_AD_data_metadata.xlsx"  # Change to your actual file path
data <- read_excel(file_path, sheet = "filtered_id_data")

# Ensure 'group' is a numeric binary variable (1 = AD, 0 = Control)
data$group <- as.numeric(data$group == "AD")  # Converts "AD" -> 1, others -> 0

# Convert categorical variables to factors
data$gender <- as.factor(data$gender)
data$race <- as.factor(data$race)
data$region <- as.factor(data$region)

# Drop ID column
data <- data %>% select(-id)

# Convert all non-numeric values (e.g., "#VALUE!") to NA
data[data == "#VALUE!"] <- NA

# Convert all numeric columns properly
numeric_columns <- setdiff(names(data), c("group", "gender", "race", "region"))
data[numeric_columns] <- lapply(data[numeric_columns], as.numeric)

# Handle missing values using mean imputation
for (col in numeric_columns) {
  data[[col]][is.na(data[[col]])] <- mean(data[[col]], na.rm = TRUE)
}

# Normalize numerical features
data[numeric_columns] <- scale(data[numeric_columns])

# Ensure 'group' is treated as numeric
print(str(data$group))  # Should show "num [1:n] 0 1 1 0 ..." (not Factor)

# Fit linear model treating AD status as continuous (Linear Probability Model)
model <- lm(group ~ ., data = data)

# Model summary to check variable significance
summary(model)

# Extract significant variables (p-value < 0.05)
significant_vars <- summary(model)$coefficients
significant_vars <- significant_vars[significant_vars[,4] < 0.05, ]
print(significant_vars)

# Sort variables by absolute coefficient size to see most predictive ones
significant_vars <- significant_vars[order(abs(significant_vars[,1]), decreasing = TRUE), ]
print(significant_vars)


#Logistric regression
set.seed(42)
trainIndex <- createDataPartition(data$group, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# Perform stepwise logistic regression
library(MASS)
model <- glm(group ~ ., data = trainData, family = binomial)
stepwise_model <- stepAIC(model, direction = "both")

# Summary of the reduced model
summary(stepwise_model)


# Fit logistic regression model
model <- glm(group ~ ., data = trainData, family = binomial)

# Model summary to check variable significance
summary(model)
