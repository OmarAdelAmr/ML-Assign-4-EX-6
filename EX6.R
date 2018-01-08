# Name: Omar Amr
# Matrikel-Nr: k11776960

library(kernlab)


directory <- '~/Desktop/ML-Assign-4-EX-6'
setwd(directory)
training_set <- read.csv('Sequences_train.csv', sep = ',', header = FALSE, stringsAsFactors = FALSE)
testing_set <- read.csv('Sequences_test_unlabeled.csv', sep = ',', header = FALSE, stringsAsFactors = FALSE)
training_data <- training_set[, 1]
training_labels <- training_set[, 2]
set.seed(123)



get_all_subsequences <- function(k, input_string)
{
  result <- substring(input_string, 1:(nchar(input_string) - k + 1), k:nchar(input_string))
  return(result)
}


count_occurances <- function(subsquence, subsequences_list)
{
  return(sum(subsequences_list == subsquence))
}


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

calculate_normalization_factor <- function(seq_length, input_string)
{
  x_subsequences <- get_all_subsequences(seq_length, input_string)
  unique_x_subsequences <- unique(x_subsequences)
  y_subsequences <- x_subsequences
  match_score <- 0
  for (sub in unique_x_subsequences)
  {
    match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
  }
  return(match_score)
}

train_spectrum_model <- function(mode = "normal", input_length = 0)
{
  if(mode == "best")
  {
    normalization_factor <- calculate_normalization_factor(3, training_set[1, 1])
    A <- matrix(nrow = matrix_size, ncol = matrix_size, byrow = TRUE)
    for (x in 1:matrix_size)
    {
      for (y in 1:matrix_size)
      {
        x_subsequences <- get_all_subsequences(3, training_set[x, 1])
        unique_x_subsequences <- unique(x_subsequences)
        y_subsequences <- get_all_subsequences(3, training_set[y, 1])
        match_score <- 0
        for (sub in unique_x_subsequences)
        {
          match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
        }
        A[x, y] <- match_score / normalization_factor
      }
    }
    A <- as.kernelMatrix(A)
    subsquence_length <<- 3
    spectrum_trained_model <<- ksvm(x=A, y=training_labels, type="C-svc", kernel='matrix', C=10, cross=10)
  } else if(mode == "search")
  {
    best_matrix <- matrix()
    min_error <- 100
    best_c <- 0
    best_sequence_length <- 0
    c_vector <- c(5, 10, 50, 100,  200)
    sequence_length_vector <- c(2, 3, 4, 5)
    
    for (c_counter in c_vector) 
    {
      print(paste0("Current C: ", c_counter))
      for (sequence_counter in sequence_length_vector)
      {
        normalization_factor <- calculate_normalization_factor(sequence_counter, training_set[1, 1])
        print(paste0("Current length: ", sequence_counter))
        A <- matrix(nrow = matrix_size, ncol = matrix_size, byrow = TRUE)
        for (x in 1:matrix_size)
        {
          for (y in 1:matrix_size)
          {
            x_subsequences <- get_all_subsequences(sequence_counter, training_set[x, 1])
            unique_x_subsequences <- unique(x_subsequences)
            y_subsequences <- get_all_subsequences(sequence_counter, training_set[y, 1])
            match_score <- 0
            for (sub in unique_x_subsequences)
            {
              match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
            }
            A[x, y] <- match_score / normalization_factor
          }
        }
        A <- as.kernelMatrix(A)
        spectrum_trained_model <<- ksvm(x=A, y=training_labels, type="C-svc", kernel='matrix', C=c_counter, cross=10)
        print(cross(spectrum_trained_model))
        if (cross(spectrum_trained_model) < min_error)
        {
          best_c <- c_counter
          best_sequence_length <- sequence_counter
          min_error <- cross(spectrum_trained_model)
          best_matrix <- A
        } 
      } 
    }
    
    print(paste0("C: ", best_c))
    print(paste0("Length: ", best_sequence_length))
    subsquence_length <<- best_sequence_length
    print(paste0("Error: ", min_error))
    spectrum_trained_model <<- ksvm(x=best_matrix, y=training_labels, type="C-svc", kernel='matrix', C=best_c, cross=10)
    
  } else if(mode == "normal")
  {
    normalization_factor <- calculate_normalization_factor(input_length, training_set[1, 1])
    A <- matrix(nrow = matrix_size, ncol = matrix_size, byrow = TRUE)
    for (x in 1:matrix_size)
    {
      for (y in 1:matrix_size)
      {
        x_subsequences <- get_all_subsequences(input_length, training_set[x, 1])
        unique_x_subsequences <- unique(x_subsequences)
        y_subsequences <- get_all_subsequences(input_length, training_set[y, 1])
        match_score <- 0
        for (sub in unique_x_subsequences)
        {
          match_score <- match_score + (count_occurances(sub, x_subsequences) * count_occurances(sub, y_subsequences))
        }
        A[x, y] <- match_score / normalization_factor
      }
    }
    A <- as.kernelMatrix(A)
    subsquence_length <<- input_length
    spectrum_trained_model <<- ksvm(x=A, y=training_labels, type="C-svc", kernel='matrix', C=10, cross=10)
  }
  
  print("Spectrum Model Training Done")
}


train_linear_model <- function()
{
  min_error <- 100
  best_c <- 0
  c_vector <- c(5, 10, 50, 100,  200)
  
  for (c_counter in c_vector) 
  {
    print(paste0("Current c: ", c_counter))
    linear_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=c_counter,kernel='vanilladot',cross=20)
    print(cross(linear_trained_model))
    if (cross(linear_trained_model) < min_error)
    {
      best_c <- c_counter
      min_error <- cross(linear_trained_model)
    }      
  }  
  print(paste0("C: ", best_c))
  print(paste0("Error: ", min_error))
  linear_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc", C=best_c, kernel='vanilladot', cross=20)
  print("Linear Model Training Done")
}


train_rbf_model <- function()
{
  min_error <- 100
  best_c <- 0
  best_sigma <- 0
  c_vector <- c(5, 10, 50, 100,  200)
  sigma_vector <- c(0.001, 0.01, 0.1, 1, 10, 100)
  
  for (c_counter in c_vector) 
  {
    for (sigma_counter in sigma_vector) 
    {
      rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=c_counter,kpar=list(sigma=sigma_counter),kernel='rbfdot',cross=50)
      print(cross(rbf_trained_model))
      if (cross(rbf_trained_model) < min_error)
      {
        best_sigma <- sigma_counter
        best_c <- c_counter
        min_error <- cross(rbf_trained_model)
      }
    }
  }
  print(paste0("Sigma: ", best_sigma))
  print(paste0("C: ", best_c))
  print(paste0("Error: ", min_error))
  rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=best_c,kpar=list(sigma=best_sigma), kernel='rbfdot', cross=50)
  print("RBF Model Training Done")
}


spectrum_model_prediction <- function(input_acid)
{
  normalization_factor <- calculate_normalization_factor(subsquence_length, training_set[1, 1])
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
    M[1, y] <- match_score / normalization_factor
  }
  M <- as.kernelMatrix(M)
  acid_prediction <- predict(spectrum_trained_model,M)
  return(acid_prediction)
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



predict_spectrum_test_sample <- function(test_data)
{
  spectrum_all_predictions_result <- vector()
  for (x in test_data) 
  {
    spectrum_all_predictions_result <- c(spectrum_all_predictions_result, spectrum_model_prediction(x))
  }
  write(spectrum_all_predictions_result, "spectrum_predictions.txt", sep="\n")
}


predict_rbf_test_sample <- function(test_data)
{
  rbf_all_predictions_result <- vector()
  for (x in test_data) 
  {
    rbf_all_predictions_result <- c(rbf_all_predictions_result, rbf_model_prediction(x))
  }
  write(rbf_all_predictions_result, "rbf_predictions.txt", sep="\n")
}


predict_linear_test_sample <- function(test_data)
{
  linear_all_predictions_result <- vector()
  for (x in test_data) 
  {
    linear_all_predictions_result <- c(linear_all_predictions_result, linear_model_prediction(x))
  }
  write(linear_all_predictions_result, "linear_predictions.txt", sep="\n")
}



subsquence_length <- 0
matrix_size <- nrow(training_set)

acid_alphabet <- vector()
acid_length <- nchar(training_data[1])
get_alphabet()
encoding_matrix <- matrix(nrow=matrix_size, ncol=acid_length*length(acid_alphabet),byrow =TRUE)
calculate_encoded_matrix()


#The whole process of testing all parameters and choosing the best to train the final model.
#################
# Train Models: # 
#################

# train_spectrum_model(mode = "search")  # search mode tests all parametrs to get the best model.
# train_spectrum_model(mode = "normal", input_length = 4)  # Set subsequence length depending on user input.
# train_linear_model()
# train_rbf_model()


# ---------------------------------------------------------- #
# Training models with paramters that acheive best results directly without search [Depending on pre-run].
################
# Best Models: # 
################

# train_spectrum_model(mode = "best") # best mode trains the model with best parameters directly without searching.
# linear_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc", C=100, kernel='vanilladot', cross=10)
rbf_trained_model <<- ksvm(encoding_matrix,training_labels,type="C-svc",C=10,kpar=list(sigma=0.001), kernel='rbfdot', cross=50)


# ---------------------------------------------------------- #
# Use trained models to predict the unlabelled test dataset. Write predictons in text file.
##################################
# Run Prediction on Test Samples #  
##################################

# predict_spectrum_test_sample(testing_set[,1])
# predict_linear_test_sample(testing_set[,1])
predict_rbf_test_sample(testing_set[,1])