library(kernlab)
library(kebabs)

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


start <- Sys.time()
train_spectrum_model()
end <- Sys.time()
print(end - start)


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
  linear_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=100, kernel='vanilladot', cross=10)
  print("Linear Model Training Done")
}


train_rbf_model <- function()
{
  rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=50,kpar=list(sigma=0.001), kernel='rbfdot', cross=10)
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
# ain_linear_model()
# train_rbf_model()
# end <- Sys.time()
# print(end - start)

# Prediction example

to_predict <- "FGEKIGLSFQLADDL"
# linear_model_prediction(to_predict)
# rbf_model_prediction(to_predict)
spectrum_model_prediction(to_predict)
