library(kernlab)

directory <- '~/Desktop/ML-Assign-4-EX-6'
setwd(directory)

training_set <- read.csv('Sequences_train.csv', sep = ',', header = FALSE, stringsAsFactors = FALSE)
testing_set <- read.csv('Sequences_test_unlabeled.csv', sep = ',', header = FALSE, stringsAsFactors = FALSE)

training_data <- training_set[, 1]
training_labels <- training_set[, 2]

subsquence_length <- 3
matrix_size <- nrow(training_set)

get_all_subsequences <- function(k, input_string)
{
  result <- substring(input_string, 1:(nchar(input_string) - k + 1), k:nchar(input_string))
  
  return(result)
}


count_occurances <- function(subsquence, subsequences_list)
{
  return(sum(subsequences_list == subsquence))
}

# ksvm
train_model <- function()
{
  
  A <- matrix(nrow = matrix_size, ncol = matrix_size, byrow = TRUE)
  
  for (x in 1:matrix_size)
  {
    for (y in 1:matrix_size)
    {
      x_subsequences <- get_all_subsequences(subsquence_length, training_set[x, 1])
      unique_x_subsequences <- unique(x_subsequences)
      y_subsequences <- get_all_subsequences(subsquence_length, training_set[y, 1])
      match_score <- 0
      for (sub in unique_x_subsequences)
      {
        match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
      }
      A[x, y] <- match_score
    }
  }
  
  A <- as(A, "KernelMatrix")
  trained_model <<- ksvm(A,training_labels,type="C-svc",C=100,scaled=c(), cross=10)
}


predict_acid <- function(input_acid)
{
  x_subsequences <- get_all_subsequences(subsquence_length, input_acid)
  M <- matrix(nrow = 1, ncol = matrix_size, byrow = TRUE)
  for (y in 1:matrix_size)
  {
    y_subsequences <- get_all_subsequences(subsquence_length, training_set[y, 1])
    match_score <- 0
    for (sub in x_subsequences)
    {
      match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
    }
    M[1, y] <- match_score
  }
  
  acid_prediction = predict(trained_model,M)
  return(acid_prediction)
}

train_model()

# -------------------------------------------------- #
# Kebabs

# library(kebabs)

# temp_train <- training_set[c(1:1800),  ]
# temp_test <- training_set[c(1801:2000),  ]
# training_labels <- temp_train[, 2]
# A <- matrix(nrow = 1800, ncol = 1800, byrow = TRUE)


# for (x in 1:1800)
# {
#   # temp_dataset <- training_set[-c(x),]
#   for (y in 1:1800)
#   {
#    
#       x_subsequences = get_all_subsequences(subsquence_length, temp_train[x, 1])
#       y_subsequences = get_all_subsequences(subsquence_length, temp_train[y, 1])
#       match_score <- 0
#       for (sub in x_subsequences)
#       {
#         match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
#       }
#       A[x, y] <- match_score
#   }
# }

# model <- kebabs:::svmd.default(A, training_labels, cost=1)

# B = matrix(nrow = 200, ncol = 1800, byrow = TRUE)
# for (x in 1:200)
# {
#   # temp_dataset <- training_set[-c(x),]
#   for (y in 1:1800)
#   {
#     x_subsequences = get_all_subsequences(subsquence_length, temp_test[x, 1])
#     y_subsequences = get_all_subsequences(subsquence_length, temp_train[y, 1])
#     match_score <- 0
#     for (sub in x_subsequences)
#     {
#       match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
#     }
#     B[x, y] <- match_score
#   }
# }
# 
# B <- as(B, "KernelMatrix")

