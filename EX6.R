library(kernlab)
library(kebabs)

set.seed(123)

start <- Sys.time()

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
train_spectrum_model <- function()
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
  spectrum_trained_model <<- ksvm(A,training_labels,type="C-svc",C=100,scaled=c(), cross=10)
  print("Spectrum Model Training Done")
}


spectrum_model_prediction <- function(input_acid)
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
  
  acid_prediction <- predict(spectrum_trained_model,M)
  return(acid_prediction)
}


# start <- Sys.time()
# train_spectrum_model()
# end <- Sys.time()
# print(end - start)


# -------------------------------------------------- #
# One hot encoding


get_alphabet <- function()
{
  for(acid in training_data)
  {
    temp <- unique(unlist(strsplit(acid, "")))
    acid_alphabet <<- unique(append(acid_alphabet, temp))
  }
  acid_alphabet <<- sort(acid_alphabet)
}


encode_sequence <- function(input_sequence)
{
  data_sample_split <- unlist(strsplit(input_sequence, ""))
  data_sample_vect <- vector()
  data_sample_encoding <- replicate(length(acid_alphabet), "0")
  for (char in data_sample_split)
  {
    temp_encoding <- replicate(length(acid_alphabet), "0")
    temp_encoding[match(char, acid_alphabet)] <- 1
    temp_encoding <- unlist(strsplit(temp_encoding, ""))
    data_sample_vect <- append(data_sample_vect, temp_encoding)
  }
  return(as.integer(data_sample_vect))
}


calculate_encoded_matrix <- function()
{
  for (data_sample_counter in 1:matrix_size)
  {
    encoding_matrix[data_sample_counter, ] <<- encode_sequence(training_data[data_sample_counter])
  }
  print("Data encoding done")
}


train_linear_model <- function()
{
  min_error <- 100

  cross_value <- 0
  c <- 0

  cross_vector <- c(10, 20, 50, 100, 200)
  c_vector <- c(5, 10, 50, 100,  200)

  for (cross_counter in cross_vector) {
    print(paste0("Current Cross: ", cross_counter))
    for (c_counter in c_vector) {
      print(paste0("Current c: ", c_counter))
        linear_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=c_counter,kernel='vanilladot',cross=cross_counter)
        print(cross(linear_trained_model))
        if (cross(linear_trained_model) < min_error)
        {
          c <- c_counter
          cross_value <- cross_counter
          min_error <- cross(linear_trained_model)
        }
    }
  }

  print(paste0("C: ", c))
  print(paste0("Cross: ", cross_value))
  
  print("Linear Model Training Done")
}


train_rbf_model <- function()
{
  min_error <- 100
  
  cross_value <- 0
  c <- 0
  sigma <- 0
  
  cross_vector <- c(10, 20, 50, 100, 200)
  c_vector <- c(5, 10, 50, 100,  200)
  sigma_vector <- c(0.001, 0.01, 0.1, 1, 10, 100)
  
  for (cross_counter in cross_vector) {
    print(paste0("Current Cross: ", cross_counter))
    for (c_counter in c_vector) {
      for (sigma_counter in sigma_vector) {
        rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=c_counter,kpar=list(sigma=sigma_counter),kernel='rbfdot',cross=cross_counter)
        print(cross(rbf_trained_model))
        if (cross(rbf_trained_model) < min_error)
        {
          sigma <- sigma_counter
          c <- c_counter
          cross_value <- cross_counter
          min_error <- cross(rbf_trained_model)
        }
      }
    }
  }
  
  # rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=50,kpar=list(sigma=0.001), kernel='rbfdot', cross=10)
  print(paste0("Sigma: ", sigma))
  print(paste0("C: ", c))
  print(paste0("Cross: ", cross_value))
  print("RBF Model Training Done")
}


linear_model_prediction <- function(input_acid)
{
  encoded <- encode_sequence(input_acid)
  prediction_mat <- matrix(nrow=1, ncol=acid_length*length(acid_alphabet),byrow =TRUE)
  prediction_mat[1, ] <- encoded
  return(predict(linear_trained_model,prediction_mat))
}

rbf_model_prediction <- function(input_acid)
{
  encoded <- encode_sequence(input_acid)
  prediction_mat <- matrix(nrow=1, ncol=acid_length*length(acid_alphabet),byrow =TRUE)
  prediction_mat[1, ] <- encoded
  return(predict(rbf_trained_model,prediction_mat))
}


acid_alphabet <- vector()
acid_length <- nchar(training_data[1])
get_alphabet()
encoding_matrix <- matrix(nrow=matrix_size, ncol=acid_length*length(acid_alphabet),byrow =TRUE)
calculate_encoded_matrix()

# start <- Sys.time()
train_linear_model()
# train_rbf_model()

# rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=50,kpar=list(sigma=0.001), kernel='rbfdot', cross=20)
# linear_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc", C=10, kernel='vanilladot', cross=20)


# TODO: To be removed
# rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=10,kpar=list(sigma=0.001), kernel='rbfdot', cross=100)

end <- Sys.time()
print(end - start)

# Prediction example

# to_predict <- "FGEKIGLSFQLADDL"
# linear_model_prediction(to_predict)
# rbf_model_prediction(to_predict)
# spectrum_model_prediction(to_predict)

# cross(rbf_trained_model)

# best rbf 0.00, 100